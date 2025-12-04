#!/usr/bin/env python3
# split_lr_presence.py
#
# Usage:
#   ./split_lr_presence.py input_wide.tsv output.vcf --original_vcf original.vcf.gz \
#       [--lr_threshold 3] [--sr_threshold 2] [--sr_only_threshold 4]
#
# Output:
#   - FILTER: TIER1 / TIER2 (else PASS)
#   - INFO: SR_ADF, SR_ADR, LR_ADF, LR_ADR, optionally ONT_ADF/ONT_ADR
#   - FORMAT + sample columns copied from original VCF when (chrom,pos,ref,alt) match

import argparse, csv, sys
import pysam
from collections import OrderedDict
from datetime import datetime
import math

def get_read_cutoffs(short_cov: int, long_cov: int,
                     error_rate: float = 0.001,
                     target_p: float = 1e-2):
    """
    Compute read count cutoffs for short, long, and combined coverage
    using Poisson approximation to the binomial distribution.

    Parameters
    ----------
    short_cov : int
        Short-read coverage
    long_cov : int
        Long-read coverage
    error_rate : float, optional
        Sequencing error rate (default = 0.001)
    target_p : float, optional
        Target probability for random error (default = 1e-2)

    Returns
    -------
    dict
        {
            'short_only_cutoff': int,
            'long_only_cutoff': int,
            'combined_sr_cutoff': int,
            'combined_lr_cutoff': int,
            'combined_total_cutoff': int
        }
    """

    def poisson_tail(lmbda, k):
        """Compute P(X >= k) for Poisson(λ)"""
        # To avoid underflow for large λ, accumulate tail from k to λ+10√λ
        # but cap at reasonable bound
        max_k = int(lmbda + 10 * math.sqrt(lmbda)) + 50
        tail = 0.0
        term = math.exp(-lmbda)  # P(X=0)
        for i in range(1, max_k + 1):
            term *= lmbda / i
            if i >= k:
                tail += term
        return tail

    def find_cutoff(cov, err_rate, target_p):
        """Find minimal count where Poisson tail < target_p."""
        lam = cov * err_rate
        for r in range(1, cov + 1):
            if poisson_tail(lam, r) < target_p:
                return r
        return cov

    # independent cutoffs
    sr_cutoff = find_cutoff(short_cov, error_rate, target_p)
    lr_cutoff = find_cutoff(long_cov, error_rate, target_p)

    # combined cutoff using total coverage
    total_cov = short_cov + long_cov
    combined_total = find_cutoff(total_cov, error_rate, target_p)

    # distribute proportionally to coverage, ensure ≥1 of each
    combined_sr = max(1, round(short_cov / total_cov * combined_total))
    combined_lr = max(1, combined_total - combined_sr)

    return {
        "short_only_cutoff": sr_cutoff,
        "long_only_cutoff": lr_cutoff,
        "combined_sr_cutoff": combined_sr,
        "combined_lr_cutoff": combined_lr,
        "combined_total_cutoff": combined_total
    }

def parse_args():
    ap = argparse.ArgumentParser(description="Tier variants by LR/SR support and preserve FORMAT/sample fields from original VCF.")
    ap.add_argument("input_tsv", help="Wide TSV from parse_minipileup.py")
    ap.add_argument("output_vcf", help="Output VCF path")
    ap.add_argument("--original_vcf", required=True, help="Original VCF with FORMAT/sample data")
    return ap.parse_args()

def safe_int(x):
    try:
        return int(x)
    except Exception:
        return 0

def load_original_records(vcf_path):
    """
    Load original VCF for header & per-record lookup.
    Return (orig_vcf, records_dict, samples)
      records_dict: (chrom,pos,ref,alt) -> pysam.VariantRecord
    """
    orig = pysam.VariantFile(vcf_path)
    records = {}
    for rec in orig.fetch():
        for alt in (rec.alts or []):
            key = (rec.chrom, rec.pos, rec.ref, alt)
            records[key] = rec
    samples = list(orig.header.samples)
    return orig, records, samples

def ensure_header_fields(hdr, ont_present):
    """Add/ensure FILTER/INFO lines exist in header."""
    # Filters
    if "TIER1" not in hdr.filters:
        hdr.add_line('##FILTER=<ID=TIER1,Description="Alt present in long-read group at or above threshold">')
    if "TIER2" not in hdr.filters:
        hdr.add_line('##FILTER=<ID=TIER2,Description="Alt present in short-read group at or above threshold; absent in long-read">')

    # INFO
    if "SR_ADF" not in hdr.info:
        hdr.add_line('##INFO=<ID=SR_ADF,Number=2,Type=Integer,Description="Short-read forward depths (REF,ALT)">')
    if "SR_ADR" not in hdr.info:
        hdr.add_line('##INFO=<ID=SR_ADR,Number=2,Type=Integer,Description="Short-read reverse depths (REF,ALT)">')
    if "LR_ADF" not in hdr.info:
        hdr.add_line('##INFO=<ID=LR_ADF,Number=2,Type=Integer,Description="Long-read forward depths (REF,ALT)">')
    if "LR_ADR" not in hdr.info:
        hdr.add_line('##INFO=<ID=LR_ADR,Number=2,Type=Integer,Description="Long-read reverse depths (REF,ALT)">')
    if ont_present:
        if "ONT_ADF" not in hdr.info:
            hdr.add_line('##INFO=<ID=ONT_ADF,Number=2,Type=Integer,Description="ONT forward depths (REF,ALT)">')
        if "ONT_ADR" not in hdr.info:
            hdr.add_line('##INFO=<ID=ONT_ADR,Number=2,Type=Integer,Description="ONT reverse depths (REF,ALT)">')

def main():
    args = parse_args()

    # Load original VCF (for header + sample data)
    orig_vcf, orig_records, samples = load_original_records(args.original_vcf)

    # Peek TSV header to detect groups (SR, LR, optional ONT)
    with open(args.input_tsv, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        header = reader.fieldnames or []
        if not header:
            sys.stderr.write("ERROR: Empty or malformed TSV header.\n")
            sys.exit(2)

        groups = OrderedDict()
        for col in header[4:]:
            if "__" in col:
                gname = col.split("__", 1)[0]
                groups.setdefault(gname, True)
        group_names = list(groups.keys())
        if len(group_names) < 2:
            sys.stderr.write("ERROR: Need at least SR and LR groups in the TSV columns.\n")
            sys.exit(2)

        sr_name, lr_name = group_names[0], group_names[1]
        ont_name = group_names[2] if len(group_names) > 2 else None

    # Build output header by copying the original
    out_header = orig_vcf.header.copy()
    out_header.add_line('##source=split_lr_presence.py')
    out_header.add_line(f'##fileDate={datetime.now().strftime("%Y%m%d")}')
    #out_header.add_line(f'##thresholds=lr:{args.lr_threshold},sr:{args.sr_threshold},sr_only:{args.sr_only_threshold}')
    ensure_header_fields(out_header, ont_present=(ont_name is not None))

    # Open output VCF with pysam (this writes a proper header)
    out_vcf = pysam.VariantFile(args.output_vcf, "w", header=out_header)

    # Process rows and emit records
    missing = 0
    with open(args.input_tsv, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")

        for row in reader:
            chrom = row["chrom"]
            pos   = int(row["pos"])
            ref   = row["ref"]
            alt   = row["alt"]

            # Counts from TSV
            sr_ref_adf = safe_int(row.get(f"{sr_name}__ref_ADF"))
            sr_ref_adr = safe_int(row.get(f"{sr_name}__ref_ADR"))
            sr_alt_adf = safe_int(row.get(f"{sr_name}__alt_ADF"))
            sr_alt_adr = safe_int(row.get(f"{sr_name}__alt_ADR"))

            lr_ref_adf = safe_int(row.get(f"{lr_name}__ref_ADF"))
            lr_ref_adr = safe_int(row.get(f"{lr_name}__ref_ADR"))
            lr_alt_adf = safe_int(row.get(f"{lr_name}__alt_ADF"))
            lr_alt_adr = safe_int(row.get(f"{lr_name}__alt_ADR"))

            ont_ref_adf = ont_ref_adr = ont_alt_adf = ont_alt_adr = 0
            if ont_name:
                ont_ref_adf = safe_int(row.get(f"{ont_name}__ref_ADF"))
                ont_ref_adr = safe_int(row.get(f"{ont_name}__ref_ADR"))
                ont_alt_adf = safe_int(row.get(f"{ont_name}__alt_ADF"))
                ont_alt_adr = safe_int(row.get(f"{ont_name}__alt_ADR"))

            sr_alt_total = sr_alt_adf + sr_alt_adr
            lr_alt_total = lr_alt_adf + lr_alt_adr

            sr_total = sr_ref_adf + sr_ref_adr + sr_alt_total 
            lr_total = lr_ref_adf + lr_ref_adr + lr_alt_total

            #print(f"row:  {row}")
            #print(f"sr_total:{sr_total}, lr_total:{lr_total}")

            if (sr_total != 0): 
                thresholds = get_read_cutoffs(sr_total,lr_total)
                # Tier classification via CLI thresholds
                if sr_alt_total >= thresholds["combined_sr_cutoff"] and lr_alt_total >= thresholds["combined_lr_cutoff"]:
                    filt = "TIER1"
                elif sr_alt_total >= thresholds["short_only_cutoff"]:
                    filt = "TIER2"
                else:
                    filt = None  
            else:
                #print("missing in minipileup")
                missing+=1
                filt = None


            # Create a new record
            rec = out_vcf.new_record(
                contig=chrom,
                start=pos - 1,                      # pysam uses 0-based start
                stop=(pos - 1) + max(1, len(ref)),  # minimal stop
                id=".",
                qual=None,
                alleles=(ref, alt)
            )

            # Set FILTER
            rec.filter.clear()
            if filt:
                rec.filter.add(filt)  # else .

            # Set INFO counts
            rec.info["SR_ADF"] = (sr_ref_adf, sr_alt_adf)
            rec.info["SR_ADR"] = (sr_ref_adr, sr_alt_adr)
            rec.info["LR_ADF"] = (lr_ref_adf, lr_alt_adf)
            rec.info["LR_ADR"] = (lr_ref_adr, lr_alt_adr)
            if ont_name:
                rec.info["ONT_ADF"] = (ont_ref_adf, ont_alt_adf)
                rec.info["ONT_ADR"] = (ont_ref_adr, ont_alt_adr)

            # Copy FORMAT/sample data from original if exact match exists
            orig_key = (chrom, pos, ref, alt)
            orig_rec = orig_records.get(orig_key)
            if orig_rec:
                # Use the same FORMAT keys in the same order as original rec
                format_keys = list(orig_rec.format.keys())
                # Ensure header already defines them (it does, since we copied header)
                for smpl in samples:
                    src_call = orig_rec.samples[smpl]
                    dst_call = rec.samples[smpl]
                    for fk in format_keys:
                        val = src_call.get(fk, None)
                        if val is not None:
                            dst_call[fk] = val
            else:
                # No matching original record: set minimal placeholder genotypes if GT exists
                if "GT" in out_header.formats:
                    for smpl in samples:
                        rec.samples[smpl]["GT"] = (None, None)

            # Write record
            out_vcf.write(rec)

    print(f"{missing} variants missing in minipileup")
    out_vcf.close()
    orig_vcf.close()

if __name__ == "__main__":
    main()

