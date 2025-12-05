#!/usr/bin/env python3
# parse_minipileup.py
#
# Usage:
#   ./parse_minipileup.py original.vcf.gz derived.vcf[.gz] output.tsv --labels SR,SR,LR,LR,ONT,ONT
#
# Requirements:
#   pip install pysam
#
# Notes:
# - original VCF: SNVs only (single REF/ALT base)
# - derived VCF: FORMAT has GT, ADF, ADR; may be multi-allelic
# - Output: one row per site; columns per label:
#     <label>__ref_ADF, <label>__ref_ADR, <label>__alt_ADF, <label>__alt_ADR
# - "ref" counts = sum of all alleles that are *not* the mosaic ALT

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

def write_caller_tags(original_vcf_path, out_path):
    targets = []
    pos_set = set()
    with open(out_path, "w") as f:
        f.write("#CHROM\tPOS\tCALLERS\n")
        with pysam.VariantFile(original_vcf_path) as vf:
            it = vf.fetch() if vf.index is not None else vf
            for rec in it:
                caller_tag = ','.join(list(rec.info["CALLERS"]))
                f.write(f"{rec.chrom}\t{rec.pos}\t{caller_tag}\n")

    
def main():
    ap = argparse.ArgumentParser(description="Extract strand-specific ref/alt counts into a wide TSV (one row per site).")
    ap.add_argument("original_vcf", help="Original SNV VCF (bgzipped/indexed recommended)")
    ap.add_argument("derived_vcf", help="Derived VCF from minipileup (with ADF/ADR). Can be bgzipped or plain; index optional.")
    ap.add_argument("out_tsv", help="Output TSV path")
    ap.add_argument("--labels", required=True,
                    help="Comma-separated labels matching samples in derived VCF, in order. "
                         "Example: SR,SR,LR,LR,ONT,ONT")
    args = ap.parse_args()

    write_caller_tags(args.original_vcf,"caller_tags.tsv")

    targets, pos_set = load_targets_from_original(args.original_vcf)
    if not targets:
        sys.stderr.write("No SNV targets found in original VCF (ensure SNV-only).\n")

    with pysam.VariantFile(args.derived_vcf) as dvf:
        sample_names = list(dvf.header.samples)

    for i in samples_names:
        lab = i.split('-')[1]
        labels.append(lab)

    if len(labels) != len(sample_names):
        sys.stderr.write(f"ERROR: {len(sample_names)} samples in VCF but {len(labels)} labels provided.\n")
        sys.exit(1)

    sample_to_label = dict(zip(sample_names, labels))
    label_set = list(dict.fromkeys(labels))  # unique labels, preserve order

    derived_by_pos = index_derived_by_position(args.derived_vcf, pos_set)

    # Prepare header
    base_cols = ["chrom", "pos", "ref", "alt"]
    per_label_cols = []
    for lab in label_set:
        per_label_cols.extend([
            f"{lab}__ref_ADF",
            f"{lab}__ref_ADR",
            f"{lab}__alt_ADF",
            f"{lab}__alt_ADR",
        ])

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(base_cols + per_label_cols) + "\n")

        for chrom, pos, ref_base, alt_base in targets:
            records = derived_by_pos.get((chrom, pos), [])

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
                    alt_idx = choose_alt_index(alleles, alt_base)

            row = [chrom, str(pos), ref_base, alt_base]

            agg = {lab: {"ref_ADF": 0, "ref_ADR": 0, "alt_ADF": 0, "alt_ADR": 0}
                   for lab in label_set}

            if chosen is not None:
                for s in sample_names:
                    call = chosen.samples.get(s, None)
                    lab = sample_to_label[s]
                    if call is None:
                        continue
                    adf = get_numberR(call, "ADF")
                    adr = get_numberR(call, "ADR")

                    # ALT counts = chosen mosaic allele
                    alt_adf = None
                    alt_adr = None
                    if alt_idx is not None:
                        try:
                            alt_adf = adf[alt_idx] if adf is not None else None
                            alt_adr = adr[alt_idx] if adr is not None else None
                        except Exception:
                            pass

                    # REF counts = sum of everything except the mosaic ALT
                    ref_adf = 0
                    ref_adr = 0
                    if adf is not None:
                        for i, v in enumerate(adf):
                            if i != alt_idx and v is not None:
                                ref_adf += v
                    if adr is not None:
                        for i, v in enumerate(adr):
                            if i != alt_idx and v is not None:
                                ref_adr += v

                    agg[lab]["ref_ADF"] += ref_adf
                    agg[lab]["ref_ADR"] += ref_adr
                    if alt_adf is not None:
                        agg[lab]["alt_ADF"] += alt_adf
                    if alt_adr is not None:
                        agg[lab]["alt_ADR"] += alt_adr

            # fill row with aggregated values
            for lab in label_set:
                row.extend([
                    str(agg[lab]["ref_ADF"]),
                    str(agg[lab]["ref_ADR"]),
                    str(agg[lab]["alt_ADF"]),
                    str(agg[lab]["alt_ADR"]),
                ])

            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()

