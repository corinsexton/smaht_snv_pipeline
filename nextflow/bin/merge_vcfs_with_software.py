#!/usr/bin/env python3
import sys
import os
import pysam
import glob
from collections import defaultdict

def detect_software(filename):
    fn = filename.lower()
    if "strelka" in fn:
        return "strelka"
    elif "mt2_to" in fn:
        return "mt2"
    else:
        return "unknown"

def merge_vcfs(sample_id, outdir="merged_vcfs"):
    tagdir = os.path.join("tagged_vcfs", sample_id)
    os.makedirs(tagdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    #vcfs = [f for f in os.listdir("/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/production_results/phasing/step5/passed") if sample_id in f and f.endswith(".vcf.gz")]
    vcfs = glob.glob(f"/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/production_results/phasing/step5/passed/*{sample_id}.phased.vcf.gz")
    print(vcfs)
    if not vcfs:
        sys.exit(f"No VCFs found for {sample_id}")

    print(f"Found {len(vcfs)} VCFs for {sample_id}:")
    for vcf in vcfs:
        print(f"  - {vcf}")

    merged_records = defaultdict(lambda: {"records": [], "softwares": set()})

    # Step 1: Read and tag each input VCF
    for vcf_path in vcfs:
        software = detect_software(vcf_path)
        if software == "unknown":
            print(f"Skipping {vcf_path} (software=unknown)")
            continue
        print(f"Reading {vcf_path} as SOFTWARE={software}")

        vcf_in = pysam.VariantFile(vcf_path)

        # Make a copy of the header and ensure SOFTWARE INFO exists
        header_copy = vcf_in.header.copy()
        if "SOFTWARE" not in header_copy.info:
            header_copy.add_line(
                '##INFO=<ID=SOFTWARE,Number=.,Type=String,Description="Source caller(s)">'
            )

        # Write tagged VCF
        tagged_path = os.path.join(
            tagdir, os.path.basename(vcf_path).replace(".vcf.gz", "_tagged.vcf.gz")
        )
        with pysam.VariantFile(tagged_path, "wz", header=header_copy) as vcf_tagged:
            for rec in vcf_in.fetch():
                new_rec = vcf_tagged.new_record(
                    contig=rec.contig,
                    start=rec.start,
                    stop=rec.stop,
                    id=rec.id,
                    qual=rec.qual,
                    alleles=rec.alleles
                )

                # Copy FILTERs
                for f in rec.filter.keys():
                    new_rec.filter.add(f)

                # Copy all INFO fields (ensure they exist in header)
                for k, v in rec.info.items():
                    if k in vcf_tagged.header.info:
                        new_rec.info[k] = v

                # Add SOFTWARE tag
                if "SOFTWARE" in new_rec.info:
                    vals = set(new_rec.info["SOFTWARE"])
                    vals.add(software)
                    new_rec.info["SOFTWARE"] = sorted(vals)
                else:
                    new_rec.info["SOFTWARE"] = [software]

                vcf_tagged.write(new_rec)

                # Collect for merged output
                key = (rec.contig, rec.pos, rec.ref, tuple(rec.alts))
                merged_records[key]["records"].append(rec)
                merged_records[key]["softwares"].add(software)

        pysam.tabix_index(tagged_path, preset="vcf", force=True)
        vcf_in.close()

    # Step 2: Prepare merged header (union of all INFO and FILTER fields)
    merged_path = os.path.join(outdir, f"{sample_id}_merged.vcf.gz")
    out_header = pysam.VariantHeader()
    out_header.add_meta("fileformat", "VCFv4.2")
    out_header.add_line(
        '##INFO=<ID=SOFTWARE,Number=.,Type=String,Description="Source caller(s)">'
    )

    info_lines = set()
    filter_lines = set()
    for vcf_path in vcfs:
        vcf_in = pysam.VariantFile(vcf_path)
        for rec in vcf_in.header.records:
            s = str(rec).strip()
            if s.startswith("##INFO=<") and "SOFTWARE" not in s:
                info_lines.add(s)
            elif s.startswith("##FILTER=<"):
                filter_lines.add(s)
        vcf_in.close()

    for line in sorted(info_lines):
        try:
            out_header.add_line(line)
        except ValueError:
            pass  # skip malformed or duplicate INFO defs

    for line in sorted(filter_lines):
        try:
            out_header.add_line(line)
        except ValueError:
            pass  # skip malformed or duplicate FILTER defs

    # Add all contigs seen
    for (chrom, _, _, _) in merged_records.keys():
        if chrom not in out_header.contigs:
            out_header.contigs.add(chrom)

    # Step 3: Write merged variants
    vcf_out = pysam.VariantFile(merged_path, "wz", header=out_header)

    print("Writing merged variants...")
    for (chrom, pos, ref, alts), data in sorted(merged_records.items()):
        rec_template = data["records"][0]
        rec = vcf_out.new_record(
            contig=chrom,
            start=pos - 1,
            stop=pos - 1 + len(ref),
            id=rec_template.id,
            qual=rec_template.qual,
            alleles=(ref,) + alts
        )

        # Merge FILTERs and INFOs
        all_filters = set()
        all_info = {}

        for r in data["records"]:
            all_filters.update(r.filter.keys())
            for k, v in r.info.items():
                if k == "SOFTWARE":
                    continue
                # Keep first value for now (simple merge)
                if k not in all_info:
                    all_info[k] = v

        for f in all_filters:
            rec.filter.add(f)

        for k, v in all_info.items():
            if k in out_header.info:
                rec.info[k] = v

        rec.info["SOFTWARE"] = sorted(data["softwares"])
        vcf_out.write(rec)

    vcf_out.close()
    pysam.tabix_index(merged_path, preset="vcf", force=True)

    print(f"Done.\n  Tagged VCFs: {tagdir}/\n  Merged result: {merged_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: merge_vcfs_with_software.py <SMHT_ID> [output_dir]")
    sample_id = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else "merged_vcfs"
    merge_vcfs(sample_id, outdir)

