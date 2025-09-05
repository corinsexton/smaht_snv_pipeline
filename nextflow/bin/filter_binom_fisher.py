#!/usr/bin/env python3
# filter_strand_and_germline.py
#
# Usage:
#   ./filter_strand_and_germline.py tiered.vcf passed.vcf \
#       [--strand_alpha 0.001] [--germline_alpha 0.001] [--min_lr_alt 2] \
#       [--min_alt_binom 2] [--summary summary.txt]
#
# Requires:
#   pip install pysam scipy
#
# Behavior:
#   - Keep only records with FILTER containing TIER1 or TIER2.
#   - Strand balance (Fisher two-sided) on chosen sample:
#       choose LR if (LR_alt_ADF + LR_alt_ADR) >= --min_lr_alt, else SR.
#       Pass if p >= --strand_alpha (i.e., NOT biased).
#   - Germline deviation (Binomial two-sided vs p=0.5):
#       run on LR if LR ALT reads >= --min_alt_binom; run on SR if SR ALT reads >= --min_alt_binom.
#       Variant passes germline if ANY available binomial test has p < --germline_alpha.
#       If no binomial tests are run (both sides below threshold), germline = FAIL.
#   - Writes only variants passing BOTH tests to output VCF.
#   - INFO annotations:
#       SB_P      : Fisher p-value (chosen sample for strand test)
#       TEST_SRC  : "LR" or "SR" used for the Fisher test
#       GLM_P_SR  : Binomial p-value on SR (if run; else omitted)
#       GLM_P_LR  : Binomial p-value on LR (if run; else omitted)
#       GLM_P     : min(GLM_P_SR, GLM_P_LR) among tests that ran

import argparse
import sys
import pysam
from datetime import datetime
from scipy.stats import fisher_exact
try:
    from scipy.stats import binomtest as _binomtest
    def binom_two_sided(k, n, p=0.5):
        return _binomtest(k, n, p, alternative="two-sided").pvalue
except Exception:
    from scipy.stats import binom_test as _binomtest
    def binom_two_sided(k, n, p=0.5):
        return _binomtest(k, n, p, alternative="two-sided")

def parse_args():
    ap = argparse.ArgumentParser(description="Filter by strand balance (Fisher) and germline deviation (binomial) using SciPy.")
    ap.add_argument("input_vcf", help="Tiered VCF from split_lr_presence.py")
    ap.add_argument("output_vcf", help="Output VCF with variants passing both tests")
    ap.add_argument("--strand_alpha", type=float, default=0.01,
                    help="Keep if Fisher p >= this (default: 0.01)")
    ap.add_argument("--germline_alpha", type=float, default=0.01,
                    help="Keep if Binomial p < this (default: 0.01)")
    ap.add_argument("--min_lr_alt", type=int, default=2,
                    help="Min LR ALT reads to use LR for the strand test (default: 2)")
    ap.add_argument("--min_alt_binom", type=int, default=1,
                    help="Min ALT reads required to run the binomial test on a sample (default: 1)")
    ap.add_argument("--summary", default=None, help="Optional summary text file")
    return ap.parse_args()

def get_pair(info, key):
    v = info.get(key, None)
    if v is None:
        return (0, 0)
    try:
        return (int(v[0]), int(v[1]))  # Number=2 -> tuple
    except Exception:
        try:
            a, b = str(v).split(",")
            return (int(a), int(b))
        except Exception:
            return (0, 0)

def main():
    args = parse_args()

    inf = pysam.VariantFile(args.input_vcf)

    # Make sure new INFO fields exist on the input header so we can set them on records.
    if "SB_P" not in inf.header.info:
        inf.header.add_line('##INFO=<ID=SB_P,Number=1,Type=Float,Description="Two-sided Fisher exact p-value for strand balance on chosen sample (LR if LR ALT >= min_lr_alt; else SR)">')
    if "GLM_P" not in inf.header.info:
        inf.header.add_line('##INFO=<ID=GLM_P,Number=1,Type=Float,Description="Minimum two-sided binomial p-value (SR/LR) among tests that ran">')
    if "GLM_P_SR" not in inf.header.info:
        inf.header.add_line('##INFO=<ID=GLM_P_SR,Number=1,Type=Float,Description="Two-sided binomial p-value vs p=0.5 on SR (if run)">')
    if "GLM_P_LR" not in inf.header.info:
        inf.header.add_line('##INFO=<ID=GLM_P_LR,Number=1,Type=Float,Description="Two-sided binomial p-value vs p=0.5 on LR (if run)">')
    if "TEST_SRC" not in inf.header.info:
        inf.header.add_line('##INFO=<ID=TEST_SRC,Number=1,Type=String,Description="Counts source used for Fisher strand test: LR or SR">')

    header = inf.header.copy()
    header.add_line('##source=filter_strand_and_germline.py')
    header.add_line(f'##fileDate={datetime.now().strftime("%Y%m%d")}')

    outf = pysam.VariantFile(args.output_vcf, "w", header=header)

    # Counters
    total = non_tier = considered = used_lr = used_sr = 0
    pass_both = fail_strand_only = fail_germ_only = fail_both = 0
    binom_tests_sr = binom_tests_lr = 0

    for rec in inf.fetch():
        total += 1
        filt_keys = set(rec.filter.keys()) if len(rec.filter.keys()) > 0 else set()
        if "TIER1" not in filt_keys and "TIER2" not in filt_keys:
            non_tier += 1
            continue
        considered += 1

        # Counts (REF,ALT) for ADF/ADR per sample
        sr_adf_ref, sr_adf_alt = get_pair(rec.info, "SR_ADF")
        sr_adr_ref, sr_adr_alt = get_pair(rec.info, "SR_ADR")
        lr_adf_ref, lr_adf_alt = get_pair(rec.info, "LR_ADF")
        lr_adr_ref, lr_adr_alt = get_pair(rec.info, "LR_ADR")

        # ----- Fisher strand test (choose LR if enough LR ALT; else SR) -----
        lr_alt_total = lr_adf_alt + lr_adr_alt
        if lr_alt_total >= args.min_lr_alt:
            src = "LR"; used_lr += 1
            ref_fwd, ref_rev = lr_adf_ref, lr_adr_ref
            alt_fwd, alt_rev = lr_adf_alt, lr_adr_alt
        else:
            src = "SR"; used_sr += 1
            ref_fwd, ref_rev = sr_adf_ref, sr_adr_ref
            alt_fwd, alt_rev = sr_adf_alt, sr_adr_alt

        try:
            _, sb_p = fisher_exact([[ref_fwd, ref_rev],[alt_fwd, alt_rev]], alternative="two-sided")
        except Exception:
            sb_p = 1.0
        strand_ok = (sb_p >= args.strand_alpha)

        # ----- Binomial (germline deviation) on BOTH samples if eligible -----
        glm_p_sr = None
        glm_p_lr = None

        # SR
        sr_alt_total = sr_adf_alt + sr_adr_alt
        sr_total_cov = sr_alt_total + (sr_adf_ref + sr_adr_ref)
        if sr_alt_total >= args.min_alt_binom and sr_total_cov > 0:
            glm_p_sr = binom_two_sided(sr_alt_total, sr_total_cov, p=0.5)
            binom_tests_sr += 1

        # LR
        lr_total_cov = lr_alt_total + (lr_adf_ref + lr_adr_ref)
        if lr_alt_total >= args.min_alt_binom and lr_total_cov > 0:
            glm_p_lr = binom_two_sided(lr_alt_total, lr_total_cov, p=0.5)
            binom_tests_lr += 1

        ran_any_binom = (glm_p_sr is not None) or (glm_p_lr is not None)
        # Pass germline if ANY available test is significant (p < alpha)
        germline_ok = False
        glm_p_min = None
        if ran_any_binom:
            candidates = []
            if glm_p_sr is not None:
                candidates.append(glm_p_sr)
            if glm_p_lr is not None:
                candidates.append(glm_p_lr)
            glm_p_min = min(candidates)
            germline_ok = any(p < args.germline_alpha for p in candidates)
        else:
            # no tests ran -> treat as not deviating from germline
            germline_ok = False

        if strand_ok and germline_ok:
            pass_both += 1
            rec.info["SB_P"] = float(sb_p)
            rec.info["TEST_SRC"] = src
            if glm_p_sr is not None:
                rec.info["GLM_P_SR"] = float(glm_p_sr)
            if glm_p_lr is not None:
                rec.info["GLM_P_LR"] = float(glm_p_lr)
            if glm_p_min is not None:
                rec.info["GLM_P"] = float(glm_p_min)
            outf.write(rec)
        else:
            if not strand_ok and not germline_ok:
                fail_both += 1
            elif not strand_ok:
                fail_strand_only += 1
            else:
                fail_germ_only += 1

    outf.close()
    inf.close()

    # Summary
    summary_lines = [
        f"Input variants:                 {total}",
        f"Dropped (non-tier):             {non_tier}",
        f"Considered (tier1/tier2):       {considered}",
        f"Used LR counts (Fisher):        {used_lr}",
        f"Used SR counts (Fisher):        {used_sr}",
        f"Binomial tests run on SR:       {binom_tests_sr}",
        f"Binomial tests run on LR:       {binom_tests_lr}",
        f"Passed both tests:              {pass_both}",
        f"Failed strand only:             {fail_strand_only}",
        f"Failed germline only:           {fail_germ_only}",
        f"Failed both:                    {fail_both}",
        f"Parameters: strand_alpha={args.strand_alpha}, germline_alpha={args.germline_alpha}, min_lr_alt={args.min_lr_alt}, min_alt_binom={args.min_alt_binom}",
    ]
    msg = "\n".join(summary_lines) + "\n"
    sys.stderr.write(msg)
    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write(msg)

if __name__ == "__main__":
    main()

