#!/usr/bin/env python3
import argparse
import os
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Sequence, Tuple, Optional

import pysam

# IUPAC mapping for the PON base
IUPAC_CONTAINS: Dict[str, set] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    # N handled specially; * handled specially
}

_tls = threading.local()

def _get_fasta(path: str) -> pysam.FastaFile:
    fa = getattr(_tls, "fasta_handle", None)
    if fa is None:
        _tls.fasta_handle = pysam.FastaFile(path)
        fa = _tls.fasta_handle
    return fa

def decider(pon_base: str, alts: Sequence[str]) -> bool:
    """
    Return True if PASS (ALT NOT present in PON), False if FAIL.
    """
    if pon_base == "*":
        return True
    if pon_base == "N" or len(pon_base) != 1:
        return False

    allowed = IUPAC_CONTAINS.get(pon_base)
    if allowed is None:
        return False

    for alt in alts or ():
        if alt is None:
            continue
        a = alt.upper()
        if len(a) == 1 and a in {"A", "C", "G", "T"} and a in allowed:
            return False
    return True

def fetch_pon_base(fasta_path: str, chrom: str, pos_1based: int) -> Optional[str]:
    try:
        fa = _get_fasta(fasta_path)
        seq = fa.fetch(chrom, pos_1based - 1, pos_1based)
        if not seq:
            return None
        return seq.upper()
    except Exception:
        return None

def choose_mode_from_ext(path: str) -> str:
    return "wz" if path.endswith(".gz") else "w"

def ensure_info_field(header: pysam.VariantHeader) -> None:
    if "BSMN_POE" not in header.info:
        header.info.add(
            "BSMN_POE",
            number=1,
            type="String",
            description="BSMN PON/PoE filter (PASS=ALT not present in PON; FAIL=ALT present/unsupported)",
        )

def clone_record_to_header(rec_in: pysam.VariantRecord, vcf_out: pysam.VariantFile) -> pysam.VariantRecord:
    """
    Create a new record bound to vcf_out.header and copy core fields, INFO, FILTER, FORMAT, samples.
    This is required if you add INFO fields to the output header.
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

    # Copy FORMAT keys/values
    for fmt_key in rec_in.format.keys():
        rec.formats[fmt_key] = rec_in.formats[fmt_key]

    # Copy per-sample fields
    for sample in rec_in.samples:
        for key, val in rec_in.samples[sample].items():
            rec.samples[sample][key] = val

    return rec

def process_record(
    idx: int,
    rec: pysam.VariantRecord,
    fasta_path: str,
    chr_regex: Optional[re.Pattern],
) -> Tuple[int, bool]:
    chrom = rec.chrom
    pos = rec.pos  # 1-based

    if chr_regex is not None and not chr_regex.match(chrom):
        return idx, False  # mimic original behavior: default FAIL

    pon_base = fetch_pon_base(fasta_path, chrom, pos)
    if pon_base is None:
        return idx, False

    alts = rec.alts or ()
    passes = decider(pon_base, alts)
    return idx, passes

def main():
    ap = argparse.ArgumentParser(description="Annotate VCF with BSMN_POE=PASS/FAIL from PON FASTA (IUPAC-aware).")
    ap.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ/BCF (indexed if compressed)")
    ap.add_argument("--fasta", required=True, help="PON FASTA (needs .fai; .gzi if bgz)")
    ap.add_argument("--out", required=True, help="Output VCF of ALL variants, annotated (.vcf or .vcf.gz)")
    ap.add_argument("--failed-out", default=None, help="Optional VCF path to also write FAIL variants (annotated).")
    ap.add_argument("--threads", type=int, default=max(os.cpu_count() or 1, 1),
                    help="Worker threads for FASTA lookups (default: CPU count)")
    ap.add_argument("--chr-regex", default=r"^chr([0-9]+|[XY])\b",
                    help="Regex contigs must match to be checked. Set '' to check all.")
    args = ap.parse_args()

    if args.failed_out and os.path.abspath(args.failed_out) == os.path.abspath(args.out):
        sys.stderr.write("[ERROR] --out and --failed-out must be different files.\n")
        sys.exit(2)

    chr_pat = re.compile(args.chr_regex) if args.chr_regex else None

    try:
        invcf = pysam.VariantFile(args.vcf)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open VCF '{args.vcf}': {e}\n")
        sys.exit(2)

    # Output header with new INFO
    out_header = invcf.header.copy()
    ensure_info_field(out_header)

    try:
        outvcf = pysam.VariantFile(args.out, choose_mode_from_ext(args.out), header=out_header)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open output '{args.out}': {e}\n")
        sys.exit(2)

    failvcf = None
    if args.failed_out:
        try:
            failvcf = pysam.VariantFile(args.failed_out, choose_mode_from_ext(args.failed_out), header=out_header)
        except Exception as e:
            sys.stderr.write(f"[ERROR] Failed to open failed-out '{args.failed_out}': {e}\n")
            sys.exit(2)

    # Warm up FASTA (fail fast if index missing)
    try:
        _get_fasta(args.fasta)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open FASTA '{args.fasta}': {e}\n")
        sys.exit(2)

    total = passed = failed = 0
    idx = 0
    pending = {}
    next_to_write = 0

    # Submit all jobs (keeps your current structure)
    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futures = []
        for rec in invcf.fetch():
            my_idx = idx
            idx += 1
            fut = ex.submit(process_record, my_idx, rec, args.fasta, chr_pat)
            futures.append((my_idx, rec, fut))

        for i, rec_in, fut in futures:
            ok = fut.result()[1]
            pending[i] = (rec_in, ok)

            while next_to_write in pending:
                r_in, ok2 = pending.pop(next_to_write)
                total += 1
                if ok2:
                    passed += 1
                    tag = "PASS"
                else:
                    failed += 1
                    tag = "FAIL"

                # Clone to output header, then set INFO
                r = clone_record_to_header(r_in, outvcf)
                r.info["BSMN_POE"] = tag

                outvcf.write(r)

                if (not ok2) and (failvcf is not None):
                    # Need a record bound to failvcf header too (same header, but safest to clone from r_in again)
                    rf = clone_record_to_header(r_in, failvcf)
                    rf.info["BSMN_POE"] = "FAIL"
                    failvcf.write(rf)

                next_to_write += 1

    outvcf.close()
    if failvcf is not None:
        failvcf.close()
    invcf.close()

    # Index outputs if .gz
    if args.out.endswith(".gz"):
        pysam.tabix_index(args.out, preset="vcf", force=True)
    if args.failed_out and args.failed_out.endswith(".gz"):
        pysam.tabix_index(args.failed_out, preset="vcf", force=True)

    sys.stderr.write(
        "[pon_annotate_vcf] Done.\n"
        f"  Input variants: {total}\n"
        f"  BSMN_POE=PASS:  {passed}\n"
        f"  BSMN_POE=FAIL:  {failed}\n"
        f"  Output VCF:     {args.out}\n"
        + (f"  FAIL VCF:       {args.failed_out}\n" if args.failed_out else "")
    )

if __name__ == "__main__":
    main()

