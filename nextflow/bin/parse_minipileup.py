#!/usr/bin/env python3
# extract_strand_counts_wide.py
#
# Usage:
#   ./extract_strand_counts_wide.py original.vcf.gz derived.vcf[.gz] output.tsv
#
# Requirements:
#   pip install pysam
#
# Notes:
# - original VCF: SNVs only (single REF/ALT base)
# - derived VCF: FORMAT has GT, ADF, ADR; may be multi-allelic; allele strings
#   can be >1 base (e.g., "AN", "GN"). For ALT mapping we use the *first base*.
# - Output: one row per site; columns per sample:
#     <sample>__ref_ADF, <sample>__ref_ADR, <sample>__alt_ADF, <sample>__alt_ADR

import sys
import argparse
import pysam
from collections import defaultdict

NA = "NA"

def first_base(allele: str):
    if allele is None:
        return None
    return allele[0] if len(allele) > 0 else None

def choose_alt_index(alleles, target_alt_base: str):
    """
    Choose ALT index (>=1) whose *first base* matches target_alt_base.
    Prefer single-base alleles. Return index or None if not found.
    """
    if target_alt_base is None or len(alleles) <= 1:
        return None
    candidates = [i for i in range(1, len(alleles)) if first_base(alleles[i]) == target_alt_base]
    if not candidates:
        return None
    singles = [i for i in candidates if len(alleles[i]) == 1]
    return singles[0] if singles else candidates[0]

def get_numberR(sample, key):
    """
    Safely get Number=R array (ADF or ADR) from a pysam Call.
    Returns list of ints or None.
    """
    val = sample.get(key, None)
    if val is None:
        return None
    try:
        return list(val)
    except Exception:
        return None

def load_targets_from_original(original_vcf_path):
    """
    Read original VCF and return:
      - targets: list of (chrom, pos, ref_base, alt_base)
      - pos_set: set of (chrom, pos)
    Assumes SNV-only in original.
    """
    targets = []
    pos_set = set()
    with pysam.VariantFile(original_vcf_path) as vf:
        it = vf.fetch() if vf.index is not None else vf
        for rec in it:
            if rec.alts is None or len(rec.alts) != 1:
                continue
            ref = rec.ref
            alt = rec.alts[0]
            if ref and alt and len(ref) == 1 and len(alt) == 1:
                targets.append((rec.chrom, rec.pos, ref, alt))
                pos_set.add((rec.chrom, rec.pos))
    return targets, pos_set

def index_derived_by_position(derived_vcf_path, pos_set):
    """
    Iterate derived VCF and keep only records at positions in pos_set.
    Return dict: (chrom, pos) -> list[pysam.VariantRecord]
    """
    by_pos = defaultdict(list)
    with pysam.VariantFile(derived_vcf_path) as vf:
        it = vf.fetch() if vf.index is not None else vf
        for rec in it:
            key = (rec.chrom, rec.pos)
            if key in pos_set:
                by_pos[key].append(rec)
    return by_pos

def main():
    ap = argparse.ArgumentParser(description="Extract strand-specific ref/alt counts into a wide TSV (one row per site).")
    ap.add_argument("original_vcf", help="Original SNV VCF (bgzipped/indexed recommended)")
    ap.add_argument("derived_vcf", help="Derived VCF from minipileup (with ADF/ADR). Can be bgzipped or plain; index optional.")
    ap.add_argument("out_tsv", help="Output TSV path")
    args = ap.parse_args()

    targets, pos_set = load_targets_from_original(args.original_vcf)
    if not targets:
        sys.stderr.write("No SNV targets found in original VCF (ensure SNV-only).\n")

    with pysam.VariantFile(args.derived_vcf) as dvf:
        sample_names = list(dvf.header.samples)

    derived_by_pos = index_derived_by_position(args.derived_vcf, pos_set)

    # Prepare header
    base_cols = ["chrom", "pos", "ref", "alt"]
    per_sample_cols = []
    for s in sample_names:
        per_sample_cols.extend([
            f"{s}__ref_ADF",
            f"{s}__ref_ADR",
            f"{s}__alt_ADF",
            f"{s}__alt_ADR",
        ])

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(base_cols + per_sample_cols) + "\n")

        for chrom, pos, ref_base, alt_base in targets:
            records = derived_by_pos.get((chrom, pos), [])

            # Pick a record at this position whose ALT set contains our original ALT (by first base).
            chosen = None
            alt_idx = None
            if records:
                for rec in records:
                    alleles = [rec.ref] + list(rec.alts or [])
                    a_idx = choose_alt_index(alleles, alt_base)
                    if a_idx is not None:
                        chosen = rec
                        alt_idx = a_idx
                        break
                if chosen is None:
                    chosen = records[0]
                    alleles = [chosen.ref] + list(chosen.alts or [])
                    alt_idx = choose_alt_index(alleles, alt_base)  # may be None
            # If we still have no record, all sample columns become NA
            row = [chrom, str(pos), ref_base, alt_base]

            if chosen is None:
                # No derived record at this site at all
                for _ in sample_names:
                    row.extend([NA, NA, NA, NA])
                out.write("\t".join(row) + "\n")
                continue

            # For REF we always use index 0; for ALT we use alt_idx (may be None)
            # Build per-sample values
            for s in sample_names:
                call = chosen.samples.get(s, None)
                if call is None:
                    row.extend([NA, NA, NA, NA])
                    continue
                adf = get_numberR(call, "ADF")
                adr = get_numberR(call, "ADR")

                def get_idx(arr, idx):
                    if arr is None or idx is None:
                        return NA
                    try:
                        return str(arr[idx])
                    except Exception:
                        return NA

                ref_adf = get_idx(adf, 0)
                ref_adr = get_idx(adr, 0)
                alt_adf = get_idx(adf, alt_idx)
                alt_adr = get_idx(adr, alt_idx)

                row.extend([ref_adf, ref_adr, alt_adf, alt_adr])

            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()

