#!/usr/bin/env python3
# split_lr_presence.py
#
# Usage:
#   ./split_lr_presence.py input_wide.tsv output.vcf [--threshold 1]
#
# Input TSV schema (from extract_strand_counts_wide.py):
#   chrom, pos, ref, alt,
#   <S1>__ref_ADF, <S1>__ref_ADR, <S1>__alt_ADF, <S1>__alt_ADR,
#   <S2>__ref_ADF, <S2>__ref_ADR, <S2>__alt_ADF, <S2>__alt_ADR, ...
#
# Semantics:
# - Sample1 (first in TSV) is SR; Sample2 (second) is LR.
# - FILTER:
#     TIER1 -> (LR_alt_ADF + LR_alt_ADR) >= threshold
#     TIER2 -> else if (SR_alt_ADF + SR_alt_ADR) >= threshold
#     .     -> otherwise (fail)
# - INFO fields (Number=2 => "REF,ALT"):
#     SR_ADF, SR_ADR, LR_ADF, LR_ADR
#   plus:
#     SR, LR (sample names)
#
# Notes:
# - Only the first two samples in the TSV are used.

import argparse
import csv
import sys
from collections import OrderedDict
from datetime import datetime

def parse_args():
    ap = argparse.ArgumentParser(description="Tier variants by alt presence in LR vs SR and emit a VCF with strand counts (REF,ALT) and sample labels in INFO.")
    ap.add_argument("input_tsv", help="Wide TSV from extract_strand_counts_wide.py")
    ap.add_argument("output_vcf", help="Output VCF path")
    return ap.parse_args()

def safe_int(x):
    try:
        return int(x)
    except Exception:
        return 0

def sanitize_info_string(s: str) -> str:
    # Keep INFO value safe (no spaces/semicolons/commas)
    return str(s).replace(";", "_").replace(" ", "_").replace(",", "_")

def main():
    args = parse_args()

    with open(args.input_tsv, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        header = reader.fieldnames or []
        if not header:
            sys.stderr.write("ERROR: Empty or malformed TSV header.\n")
            sys.exit(2)

        # Detect samples from "<sample>__..." columns; use the first two as SR, LR.
        samples = OrderedDict()
        for col in header[4:]:
            if "__" in col:
                sname = col.split("__", 1)[0]
                samples.setdefault(sname, True)

        sample_names = list(samples.keys())
        if len(sample_names) < 2:
            sys.stderr.write("ERROR: Need at least two samples in TSV header.\n")
            sys.exit(2)

        sr_name, lr_name = sample_names[0], sample_names[1]

        # Required columns
        needed = [
            "chrom", "pos", "ref", "alt",
            f"{sr_name}__ref_ADF", f"{sr_name}__ref_ADR", f"{sr_name}__alt_ADF", f"{sr_name}__alt_ADR",
            f"{lr_name}__ref_ADF", f"{lr_name}__ref_ADR", f"{lr_name}__alt_ADF", f"{lr_name}__alt_ADR",
        ]
        missing = [c for c in needed if c not in header]
        if missing:
            sys.stderr.write("ERROR: Missing expected columns: " + ", ".join(missing) + "\n")
            sys.exit(2)

        with open(args.output_vcf, "w") as fout:
            # VCF header
            fout.write("##fileformat=VCFv4.2\n")
            fout.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            fout.write("##source=split_lr_presence.py\n")
            fout.write('##FILTER=<ID=TIER1,Description="Alt present in second (long-read) sample at or above threshold">\n')
            fout.write('##FILTER=<ID=TIER2,Description="Alt present in first (short-read) sample at or above threshold; absent in second">\n')
            fout.write('##INFO=<ID=SR_ADF,Number=2,Type=Integer,Description="Short-read sample forward depths (ADF) as REF,ALT">\n')
            fout.write('##INFO=<ID=SR_ADR,Number=2,Type=Integer,Description="Short-read sample reverse depths (ADR) as REF,ALT">\n')
            fout.write('##INFO=<ID=LR_ADF,Number=2,Type=Integer,Description="Long-read sample forward depths (ADF) as REF,ALT">\n')
            fout.write('##INFO=<ID=LR_ADR,Number=2,Type=Integer,Description="Long-read sample reverse depths (ADR) as REF,ALT">\n')
            fout.write('##INFO=<ID=SR,Number=1,Type=String,Description="Short-read sample name">\n')
            fout.write('##INFO=<ID=LR,Number=1,Type=String,Description="Long-read sample name">\n')
            fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            for row in reader:
                chrom = row["chrom"]
                pos = row["pos"]
                ref = row["ref"]
                alt = row["alt"]

                # SR counts
                sr_ref_adf = safe_int(row.get(f"{sr_name}__ref_ADF"))
                sr_ref_adr = safe_int(row.get(f"{sr_name}__ref_ADR"))
                sr_alt_adf = safe_int(row.get(f"{sr_name}__alt_ADF"))
                sr_alt_adr = safe_int(row.get(f"{sr_name}__alt_ADR"))

                # LR counts
                lr_ref_adf = safe_int(row.get(f"{lr_name}__ref_ADF"))
                lr_ref_adr = safe_int(row.get(f"{lr_name}__ref_ADR"))
                lr_alt_adf = safe_int(row.get(f"{lr_name}__alt_ADF"))
                lr_alt_adr = safe_int(row.get(f"{lr_name}__alt_ADR"))

                # Tiers use ALT totals only
                sr_alt_total = sr_alt_adf + sr_alt_adr
                lr_alt_total = lr_alt_adf + lr_alt_adr

                lr_threshold = 1
                sr_threshold = 2

                if lr_alt_total >= lr_threshold:
                    filt = "TIER1"
                elif sr_alt_total >= sr_threshold:
                    filt = "TIER2"
                else:
                    filt = "."

                info = (
                    f"SR_ADF={sr_ref_adf},{sr_alt_adf};"
                    f"SR_ADR={sr_ref_adr},{sr_alt_adr};"
                    f"LR_ADF={lr_ref_adf},{lr_alt_adf};"
                    f"LR_ADR={lr_ref_adr},{lr_alt_adr};"
                    f"SR={sanitize_info_string(sr_name)};LR={sanitize_info_string(lr_name)}"
                )

                # ID and QUAL set to '.'
                fout.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}\n")

if __name__ == "__main__":
    main()

