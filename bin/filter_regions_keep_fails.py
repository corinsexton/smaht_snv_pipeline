#!/usr/bin/env python3
import argparse
import bisect
import gzip
import os
from typing import Dict, List, Optional, Tuple

import pysam


# ----------------------------
# Interval index (0-based, half-open [start, end))
# ----------------------------

class IntervalIndex:
    def __init__(self) -> None:
        self._starts: Dict[str, List[int]] = {}
        self._ends: Dict[str, List[int]] = {}

    def add(self, chrom: str, start0: int, end0: int) -> None:
        if end0 <= start0:
            return
        self._starts.setdefault(chrom, []).append(start0)
        self._ends.setdefault(chrom, []).append(end0)

    def finalize_merge(self) -> None:
        """Sort and merge overlaps per chrom (in-place)."""
        for chrom in list(self._starts.keys()):
            starts = self._starts[chrom]
            ends = self._ends[chrom]
            order = sorted(range(len(starts)), key=lambda i: starts[i])
            starts = [starts[i] for i in order]
            ends = [ends[i] for i in order]

            merged_s: List[int] = []
            merged_e: List[int] = []
            for s, e in zip(starts, ends):
                if not merged_s or s > merged_e[-1]:
                    merged_s.append(s)
                    merged_e.append(e)
                else:
                    if e > merged_e[-1]:
                        merged_e[-1] = e

            self._starts[chrom] = merged_s
            self._ends[chrom] = merged_e

    def overlaps(self, chrom: str, start0: int, end0: int) -> bool:
        s_list = self._starts.get(chrom)
        if not s_list:
            return False
        e_list = self._ends[chrom]

        # Find the rightmost interval with start < end0
        i = bisect.bisect_left(s_list, end0) - 1
        if i < 0:
            return False

        # Because merged & sorted, checking this one and maybe earlier is enough.
        # Overlap condition (half-open): interval_end > start0
        while i >= 0 and s_list[i] < end0:
            if e_list[i] > start0:
                return True
            i -= 1
        return False


# ----------------------------
# Region sources
# ----------------------------

class RegionSource:
    label: str
    def overlaps(self, chrom: str, start0: int, end0: int) -> bool:
        raise NotImplementedError


class TabixBedRegionSource(RegionSource):
    def __init__(self, path: str, label: str):
        self.path = path
        self.label = label
        self.tbx = pysam.TabixFile(path)

    def overlaps(self, chrom: str, start0: int, end0: int) -> bool:
        try:
            it = self.tbx.fetch(chrom, start0, end0)
            for _ in it:
                return True
            return False
        except ValueError:
            return False


class InMemoryBedRegionSource(RegionSource):
    def __init__(self, bed_path: str, label: str):
        self.path = bed_path
        self.label = label
        self.idx = IntervalIndex()
        self._load_bed(bed_path)
        self.idx.finalize_merge()

    def _load_bed(self, bed_path: str) -> None:
        opener = open
        if bed_path.endswith(".gz"):
            opener = gzip.open
        with opener(bed_path, "rt") as f:
            for line in f:
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    s = int(parts[1])
                    e = int(parts[2])
                except ValueError:
                    continue
                self.idx.add(chrom, s, e)

    def overlaps(self, chrom: str, start0: int, end0: int) -> bool:
        return self.idx.overlaps(chrom, start0, end0)


class OneKgIndelsSlopRegionSource(RegionSource):
    """
    Build merged intervals from a 1KG indels VCF by expanding each indel by +/- slop bp.
    Equivalent to:
      bcftools query -f'%CHROM\t%POS0\t%END\n' kg_indels.vcf.gz | slop +/-5 | bedtools merge
    """
    def __init__(self, kg_indels_vcf: str, label: str = "1kg_indels", slop: int = 5):
        self.path = kg_indels_vcf
        self.label = label
        self.slop = slop
        self.idx = IntervalIndex()
        self._load_kg_indels_vcf(kg_indels_vcf, slop)
        self.idx.finalize_merge()

    def _load_kg_indels_vcf(self, vcf_path: str, slop: int) -> None:
        vcf = pysam.VariantFile(vcf_path)
        for rec in vcf.fetch():
            # pysam: rec.pos is 1-based; rec.stop is 1-based inclusive end in many contexts.
            # For interval work: [pos-1, stop) is a safe half-open span for the allele.
            start0 = rec.pos - 1
            end0 = rec.stop

            s = start0 - slop
            if s < 0:
                s = 0
            e = end0 + slop
            self.idx.add(rec.chrom, s, e)
        vcf.close()

    def overlaps(self, chrom: str, start0: int, end0: int) -> bool:
        return self.idx.overlaps(chrom, start0, end0)


def open_bed_region_source(path: Optional[str], label: str) -> Optional[RegionSource]:
    if not path:
        return None
    # Prefer tabix if present (bgz bed + .tbi)
    if path.endswith(".gz") and os.path.exists(path + ".tbi"):
        return TabixBedRegionSource(path, label)
    return InMemoryBedRegionSource(path, label)


# ----------------------------
# VCF annotation
# ----------------------------

def ensure_info_fields(header: pysam.VariantHeader) -> None:
    if "REGIONS" not in header.info:
        header.info.add(
            "REGIONS",
            number=1,
            type="String",
            description="Region filter annotation (PASS=not in excluded regions, FAIL=overlaps excluded regions)",
        )
    if "REGIONS_FAILED" not in header.info:
        header.info.add(
            "REGIONS_FAILED",
            number=".",
            type="String",
            description="If REGIONS=FAIL, list of region sets overlapped (SegDup,Centr,Simple_Repeat,1kg_indels)",
        )


def clone_record_to_header(rec_in: pysam.VariantRecord, vcf_out: pysam.VariantFile) -> pysam.VariantRecord:
    """
    Create a new record bound to vcf_out.header and copy core fields, INFO, FILTER, FORMAT, samples.
    This avoids 'unknown INFO' when adding new tags present only in output header.
    """
    rec = vcf_out.new_record(
        contig=rec_in.contig,
        start=rec_in.start,
        stop=rec_in.stop,
        id=rec_in.id,
        alleles=rec_in.alleles,
        qual=rec_in.qual,
        filter=list(rec_in.filter.keys()) if rec_in.filter is not None else None,
        info=dict(rec_in.info),
    )

    # Copy FORMAT values
    for fmt_key in rec_in.format.keys():
        rec.formats[fmt_key] = rec_in.formats[fmt_key]

    # Copy per-sample fields
    for sample in rec_in.samples:
        for key, val in rec_in.samples[sample].items():
            rec.samples[sample][key] = val

    return rec



def annotate_regions(
    input_vcf: str,
    output_vcf: str,
    output_vcf_failed: str,
    segdup_bed: Optional[str],
    centromere_bed: Optional[str],
    simple_repeat_bed: Optional[str],
    kg_indels_vcf: Optional[str],
    kg_slop: int = 5,
) -> None:
    vcf_in = pysam.VariantFile(input_vcf)

    header = vcf_in.header.copy()
    ensure_info_fields(header)

    mode = "wz" if output_vcf.endswith(".gz") else "w"
    vcf_out = pysam.VariantFile(output_vcf, mode, header=header)
    vcf_out_failed = pysam.VariantFile(output_vcf_failed, mode, header=header)

    sources: List[RegionSource] = []

    s = open_bed_region_source(segdup_bed, "SegDup")
    if s:
        sources.append(s)
    c = open_bed_region_source(centromere_bed, "Centr")
    if c:
        sources.append(c)
    r = open_bed_region_source(simple_repeat_bed, "Simple_Repeat")
    if r:
        sources.append(r)

    if kg_indels_vcf:
        sources.append(OneKgIndelsSlopRegionSource(kg_indels_vcf, "1kg_indels", slop=kg_slop))

    for rec_in in vcf_in.fetch():
        rec = clone_record_to_header(rec_in, vcf_out)

        # Variant span in 0-based half-open coordinates
        start0 = rec.pos - 1
        end0 = rec.stop

        failed: List[str] = []
        for src in sources:
            if src.overlaps(rec.chrom, start0, end0):
                failed.append(src.label)

        if failed:
            rec.info["REGIONS"] = "FAIL"
            rec.info["REGIONS_FAILED"] = failed  # pysam writes Number=. as comma-separated
            vcf_out_failed.write(rec)
        else:
            rec.info["REGIONS"] = "PASS"
            if "REGIONS_FAILED" in rec.info:
                del rec.info["REGIONS_FAILED"]
            vcf_out.write(rec)


    vcf_in.close()
    vcf_out.close()
    vcf_out_failed.close()

    if output_vcf.endswith(".gz"):
        pysam.tabix_index(output_vcf, preset="vcf", force=True)
        pysam.tabix_index(output_vcf_failed, preset="vcf", force=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate VCF records with REGIONS=PASS/FAIL and REGIONS_FAILED based on overlaps with region sets."
    )
    parser.add_argument("input_vcf", help="Input VCF (bgzipped or not)")
    parser.add_argument("output_vcf", help="Output VCF (bgzipped if ends with .gz)")
    parser.add_argument("output_vcf_failed", help="Output VCF failed vars (bgzipped if ends with .gz)")

    parser.add_argument("--segdup", default=None, help="UCSC SegDup BED (bgz + .tbi preferred)")
    parser.add_argument("--centromere", default=None, help="Centromere BED (bgz + .tbi preferred)")
    parser.add_argument("--simple-repeat", default=None, help="Simple repeat BED (bgz + .tbi preferred)")

    parser.add_argument("--kg-indels", default=None, help="1KG indels VCF (bgz + .tbi preferred)")
    parser.add_argument("--kg-slop", type=int, default=5, help="Slop bp for 1KG indels expansion [default=5]")

    args = parser.parse_args()

    annotate_regions(
        args.input_vcf,
        args.output_vcf,
        args.output_vcf_failed,
        segdup_bed=args.segdup,
        centromere_bed=args.centromere,
        simple_repeat_bed=args.simple_repeat,
        kg_indels_vcf=args.kg_indels,
        kg_slop=args.kg_slop,
    )

