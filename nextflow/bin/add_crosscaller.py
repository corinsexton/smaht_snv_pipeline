#!/usr/bin/env python3

import argparse
import pysam
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add CrossCaller flag when INFO/CALLERS has more than one value."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input VCF (bgzipped or uncompressed)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output VCF (bgzipped recommended)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Open input VCF
    try:
        vcf_in = pysam.VariantFile(args.input)
    except Exception as e:
        sys.exit(f"ERROR: Failed to open input VCF: {e}")

    # Copy header
    header = vcf_in.header.copy()

    # Open output
    mode = "w"
    if not args.output.endswith(".gz"):
        sys.stderr.write(
            "WARNING: Output VCF is not bgzipped. Consider using .vcf.gz for indexing.\n"
        )

    vcf_out = pysam.VariantFile(args.output, mode, header=header)

    # Process records
    for rec in vcf_in:
        if len(rec.info.get("CALLERS")) > 1:
            rec.info["CrossCaller"] = True  # Add flag

        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()

