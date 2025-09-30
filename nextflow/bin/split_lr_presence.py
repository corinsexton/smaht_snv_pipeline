#!/usr/bin/env python3
# split_lr_presence.py
#
# Usage:
#   ./split_lr_presence.py input_wide.tsv output.vcf
#
# Input TSV schema (from extract_strand_counts_wide.py):
#   chrom, pos, ref, alt,
#   <SR>__ref_ADF, <SR>__ref_ADR, <SR>__alt_ADF, <SR>__alt_ADR,
#   <LR>__ref_ADF, <LR>__ref_ADR, <LR>__alt_ADF, <LR>__alt_ADR,
#   [optional ONT group columns ...]
#
# Semantics:
# - First label in TSV → SR
# - Second label in TSV → LR
# - Third label in TSV (if present) → ONT
#
# FILTER:
#   TIER1 -> (LR_alt_ADF + LR_alt_ADR) >= 2
#   TIER2 -> else if (SR_alt_ADF + SR_alt_ADR) >= 2
#   .     -> otherwise
#
# INFO fields (Number=2 = "REF,ALT"):
#   SR_ADF, SR_ADR, LR_ADF, LR_ADR
#   ONT_ADF, ONT_ADR (if ONT present)
#   SR, LR, ONT (sample/group names)

import argparse
import csv
import sys
from collections import OrderedDict
from datetime import datetime


def parse_args():
    ap = argparse.ArgumentParser(description="Tier variants by alt presence in LR vs SR and emit a VCF with strand counts (REF,ALT) and group labels in INFO.")
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

        # Detect groups from "<group>__..." columns
        groups = OrderedDict()
        for col in header[4:]:
            if "__" in col:
                gname = col.split("__", 1)[0]
                groups.setdefault(gname, True)

        group_names = list(groups.keys())
        if len(group_names) < 2:
            sys.stderr.write("ERROR: Need at least SR and LR groups in TSV header.\n")
            sys.exit(2)

        sr_name, lr_name = group_names[0], group_names[1]
        ont_name = group_names[2] if len(group_names) > 2 else None

        # Required columns
        needed = [
            "chrom", "pos", "ref", "alt",
            f"{sr_name}__ref_ADF", f"{sr_name}__ref_ADR", f"{sr_name}__alt_ADF", f"{sr_name}__alt_ADR",
            f"{lr_name}__ref_ADF", f"{lr_name}__ref_ADR", f"{lr_name}__alt_ADF", f"{lr_name}__alt_ADR",
        ]
        if ont_name:
            needed.extend([
                f"{ont_name}__ref_ADF", f"{ont_name}__ref_ADR", f"{ont_name}__alt_ADF", f"{ont_name}__alt_ADR"
            ])

        missing = [c for c in needed if c not in header]
        if missing:
            sys.stderr.write("ERROR: Missing expected columns: " + ", ".join(missing) + "\n")
            sys.exit(2)

        with open(args.output_vcf, "w") as fout:
            # VCF header
            fout.write("##fileformat=VCFv4.2\n")
            fout.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            fout.write("##source=split_lr_presence.py\n")
            fout.write('##FILTER=<ID=TIER1,Description="Alt present in long-read group at or above threshold">\n')
            fout.write('##FILTER=<ID=TIER2,Description="Alt present in short-read group at or above threshold; absent in long-read">\n')
            fout.write('##INFO=<ID=SR_ADF,Number=2,Type=Integer,Description="Short-read forward depths (REF,ALT)">\n')
            fout.write('##INFO=<ID=SR_ADR,Number=2,Type=Integer,Description="Short-read reverse depths (REF,ALT)">\n')
            fout.write('##INFO=<ID=LR_ADF,Number=2,Type=Integer,Description="Long-read forward depths (REF,ALT)">\n')
            fout.write('##INFO=<ID=LR_ADR,Number=2,Type=Integer,Description="Long-read reverse depths (REF,ALT)">\n')
            if ont_name:
                fout.write('##INFO=<ID=ONT_ADF,Number=2,Type=Integer,Description="ONT forward depths (REF,ALT)">\n')
                fout.write('##INFO=<ID=ONT_ADR,Number=2,Type=Integer,Description="ONT reverse depths (REF,ALT)">\n')
            fout.write('##INFO=<ID=SR,Number=1,Type=String,Description="Short-read group label">\n')
            fout.write('##INFO=<ID=LR,Number=1,Type=String,Description="Long-read group label">\n')
            if ont_name:
                fout.write('##INFO=<ID=ONT,Number=1,Type=String,Description="ONT group label">\n')
            fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            for row in reader:
                chrom, pos, ref, alt = row["chrom"], row["pos"], row["ref"], row["alt"]

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

                # ONT counts (optional)
                ont_ref_adf = ont_ref_adr = ont_alt_adf = ont_alt_adr = 0
                if ont_name:
                    ont_ref_adf = safe_int(row.get(f"{ont_name}__ref_ADF"))
                    ont_ref_adr = safe_int(row.get(f"{ont_name}__ref_ADR"))
                    ont_alt_adf = safe_int(row.get(f"{ont_name}__alt_ADF"))
                    ont_alt_adr = safe_int(row.get(f"{ont_name}__alt_ADR"))

                # Apply thresholds
                sr_alt_total = sr_alt_adf + sr_alt_adr
                lr_alt_total = lr_alt_adf + lr_alt_adr


                # FROM CHATGPT
                #   Recommendations
                #   When you have PacBio 600×:
                #       Minimum: Illumina ≥2 supporting reads and HiFi ≥3 supporting reads.
                #           Expected random co-validation ≈ 0.053 per 10k candidates → already very low.
                #       Conservative: Illumina ≥2 and HiFi ≥4 → ~0.0026 per 10k.
                #           If you expect very large candidate sets (≥100k), bump either Illumina to ≥3 or HiFi to ≥4 to keep randoms ≪1.
                #   Without any PacBio support (Illumina-only):
                #       Practical default: Illumina ≥4 reads (expected random ≈ 0.038 per 10k).
                #       Large candidate sets (≥100k) or clinical-stringent: Illumina ≥5 reads (≈ 0.007 per 100k).
                #       Avoid relying on ≥3 reads alone unless you’re filtering down to ≪10k candidates; ≥3 gives ≈ 1.5 per 10k randoms.
                lr_threshold = 3 
                sr_threshold = 2 

                sr_only_threshold = 4

                if sr_alt_total >=  sr_threshold and lr_alt_total >= lr_threshold:
                    filt = "TIER1"
                elif sr_alt_total >= sr_only_threshold:
                    filt = "TIER2"
                else:
                    filt = "."

                # Build INFO string
                info_parts = [
                    f"SR_ADF={sr_ref_adf},{sr_alt_adf}",
                    f"SR_ADR={sr_ref_adr},{sr_alt_adr}",
                    f"LR_ADF={lr_ref_adf},{lr_alt_adf}",
                    f"LR_ADR={lr_ref_adr},{lr_alt_adr}",
                    f"SR={sanitize_info_string(sr_name)}",
                    f"LR={sanitize_info_string(lr_name)}",
                ]
                if ont_name:
                    info_parts.extend([
                        f"ONT_ADF={ont_ref_adf},{ont_alt_adf}",
                        f"ONT_ADR={ont_ref_adr},{ont_alt_adr}",
                        f"ONT={sanitize_info_string(ont_name)}",
                    ])

                info = ";".join(info_parts)
                fout.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}\n")


if __name__ == "__main__":
    main()

