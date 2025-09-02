#!/usr/bin/env python3
# vcf_depth_filter.py
# Usage: ./vcf_depth_filter.py input.vcf[.gz] sample.bam output.vcf[.gz] [--threshold 300]

import argparse
import sys
import pysam
import os

def main():
    p = argparse.ArgumentParser(description="Filter VCF by mean depth across REF span using pysam.")
    p.add_argument("vcf", help="Input VCF (optionally bgzipped)")
    p.add_argument("bam", help="Input BAM (indexed)")
    p.add_argument("out", help="Output VCF (use .gz to bgzip)")
    p.add_argument("--threshold", type=float, default=300.0,
                   help="Keep variants with mean depth < THRESHOLD (default: 300)")
    p.add_argument("--baseq", type=int, default=0,
                   help="Base quality threshold for counting coverage (default: 0)")
    args = p.parse_args()

    # Open BAM
    try:
        bam = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        sys.stderr.write(f"ERROR opening BAM: {e}\n")
        sys.exit(1)

    # Open VCF in/out
    try:
        vcf_in = pysam.VariantFile(args.vcf)
    except Exception as e:
        sys.stderr.write(f"ERROR opening VCF: {e}\n")
        sys.exit(1)

    header = vcf_in.header.copy()
    # Add INFO for mean depth if not present
    try:
        header.add_line('##INFO=<ID=DP_AVG,Number=1,Type=Float,Description="Mean depth across REF span computed from BAM">')
    except ValueError:
        # Already present; fine to proceed
        pass

    out_mode = "wz" if args.out.endswith(".gz") else "w"
    try:
        vcf_out = pysam.VariantFile(args.out, out_mode, header=header)
    except Exception as e:
        sys.stderr.write(f"ERROR creating output VCF: {e}\n")
        sys.exit(1)

    kept = 0
    total = 0

    for rec in vcf_in:
        total += 1
        chrom = rec.chrom
        pos1  = rec.pos                 # 1-based
        ref   = rec.ref
        ref_len = len(ref)

        # Convert to 0-based, end-exclusive for pysam
        start0 = pos1 - 1
        end0   = start0 + ref_len

        try:
            # Returns 4 arrays (A,C,G,T) of per-base counts; sum = depth
            A,C,G,T = bam.count_coverage(chrom, start0, end0,
                                         quality_threshold=args.baseq,
                                         read_callback="all")
        except ValueError as e:
            # Typically contig not found or region out of bounds
            sys.stderr.write(f"Warning: {chrom}:{pos1} skipped ({e})\n")
            continue

        npos = end0 - start0
        if npos <= 0:
            sys.stderr.write(f"Warning: {chrom}:{pos1} invalid span; skipping\n")
            continue

        depths = [A[i] + C[i] + G[i] + T[i] for i in range(npos)]
        if not depths:
            continue

        mean_depth = sum(depths) / float(len(depths))

        # Store the computed depth
        rec.info["DP_AVG"] = round(mean_depth, 2)

        if mean_depth < args.threshold:
            vcf_out.write(rec)
            kept += 1

    vcf_out.close()
    bam.close()
    vcf_in.close()

    sys.stderr.write(f"Done. Kept {kept} / {total} variants with DP_AVG < {args.threshold}.\n")

if __name__ == "__main__":
    main()

