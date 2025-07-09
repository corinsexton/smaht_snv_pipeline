#!/usr/bin/env python

import pysam
import argparse

def get_variant_key(record):
    return f"{record.chrom}:{record.pos}:{record.ref}:{','.join(record.alts)}"

def filter_vcf(input_vcf, error_vcf, output_vcf):
    vcf_in = pysam.VariantFile(input_vcf)
    panel_in = pysam.VariantFile(error_vcf)
    vcf_out = pysam.VariantFile(output_vcf, "wz", header=vcf_in.header)

    panel_keys = {
        get_variant_key(rec)
        for rec in panel_in
    }

    for rec in vcf_in:
        if get_variant_key(rec) not in panel_keys:
            vcf_out.write(rec)

    vcf_in.close()
    panel_in.close()
    vcf_out.close()

    pysam.tabix_index(output_vcf, preset="vcf", force=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove variants matching a known error panel VCF")
    parser.add_argument("input_vcf", help="Input VCF to filter (bgzipped)")
    parser.add_argument("error_vcf", help="Error panel VCF (bgzipped)")
    parser.add_argument("output_vcf", help="Output filtered VCF (bgzipped)")

    args = parser.parse_args()
    filter_vcf(args.input_vcf, args.error_vcf, args.output_vcf)

