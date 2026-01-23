#!/usr/bin/env python3
import argparse
import pysam
from collections import defaultdict

############################################################
# Parse arguments
############################################################
parser = argparse.ArgumentParser(description="Annotate VCF using Minipileup results")
parser.add_argument("--tissue", required=True, help="Current tissue (e.g. 3A)")
parser.add_argument("--orig_vcf", required=True, help="Original input VCF (bgzipped)")
parser.add_argument("--mp_vcf", required=True, help="Minipileup VCF (bgzipped)")
parser.add_argument("--out", required=True, help="Output annotated VCF file")
args = parser.parse_args()

current_tissue = args.tissue

############################################################
# 1. Load mp_vcf and compute sample → tissue VAFs
############################################################
mp = pysam.VariantFile(args.mp_vcf)

sample2tissue = {s: s.split("-")[-1] for s in mp.header.samples}

# Store VAFs for each sample, then aggregate per tissue
tissue_vafs = defaultdict(list)  # key = (chrom,pos,tissue) → list of vafs

for rec in mp:
    chrom = rec.chrom
    pos = rec.pos

    for sample in rec.samples:
        call = rec.samples[sample]

        if "ADF" not in call or "ADR" not in call:
            continue

        adf = call["ADF"]
        adr = call["ADR"]

        ref_count = adf[0] + adr[0]
        alt_count = adf[1] + adr[1]
        depth = ref_count + alt_count
        if depth == 0:
            continue

        vaf = alt_count / depth
        tissue = sample2tissue[sample]

        tissue_vafs[(chrom, pos, tissue)].append(vaf)

############################################################
# 2. Aggregate per tissue and REMOVE 0-VAF tissues
############################################################
aggregated_vaf = defaultdict(dict)   # (chrom,pos) → {tissue: vaf}

for (chrom, pos, tissue), vafs in tissue_vafs.items():
    mean_vaf = sum(vafs) / len(vafs)

    #  Skip tissues with zero VAF
    if mean_vaf > 0:
        aggregated_vaf[(chrom, pos)][tissue] = mean_vaf

############################################################
# 3. Build TISSUE_SR_VAFS and CrossTissue annotations
############################################################
summary = {}           # (chrom,pos) → "3A,0.05;3I,0.2"
crosstissue_flag = {}  # variants shared >1 tissue

for key, tissue_dict in aggregated_vaf.items():
    chrom, pos = key

    # skip if no tissues remain after filtering
    if len(tissue_dict) == 0:
        continue

    # Build tissue strings
    parts = [f"{t},{vaf:.6f}" for t, vaf in tissue_dict.items()]
    summary[key] = "|".join(parts)

    # CrossTissue if more than one tissue has VAF
    if len(tissue_dict) > 1:
        crosstissue_flag[key] = True

############################################################
# 4. Annotate the original VCF
############################################################
orig = pysam.VariantFile(args.orig_vcf)

orig.header.info.add("SR_VAF", number=1, type="Float",
                     description="VAF for short read in current tissue")
orig.header.info.add("PB_VAF", number=1, type="Float",
                     description="VAF for PacBio in current donor pooled tissues")
orig.header.info.add("ONT_VAF", number=1, type="Float",
                     description="VAF for ONT in current donor pooled tissues")
orig.header.info.add("TISSUE_SR_VAFS", number=".", type="String",
                     description="VAFs for all tissues with nonzero VAF")
orig.header.info.add("CrossTissue", number=0, type="Flag",
                     description="Variant has VAF > 0 in another tissue")

out = pysam.VariantFile(args.out, "w", header=orig.header)

############################################################
# 5. Write output with annotations
############################################################
for rec in orig:
    key = (rec.chrom, rec.pos)

    # Add current tissue VAF only if >0
    if current_tissue in aggregated_vaf.get(key, {}):
        rec.info["SR_VAF"] = float(aggregated_vaf[key][current_tissue])

    # Add TISSUE_SR_VAFS if any
    if key in summary:
        rec.info["TISSUE_SR_VAFS"] = summary[key]

    # Add CrossTissue flag
    if key in crosstissue_flag:
        rec.info["CrossTissue"] = True

    ############################################################
    # NEW: Calculate PB_VAF (using LR_ADF/LR_ADR)
    ############################################################
    try:
        lr_adf = rec.info.get("PB_ADF")
        lr_adr = rec.info.get("PB_ADR")
        if lr_adf and lr_adr:
            ref = lr_adf[0] + lr_adr[0]
            alt = lr_adf[1] + lr_adr[1]
            if ref + alt > 0:
                rec.info["PB_VAF"] = float(alt / (ref + alt))
    except Exception:
        pass

    ############################################################
    # NEW: Calculate ONT_VAF (using ONT_ADF/ONT_ADR)
    ############################################################
    try:
        ont_adf = rec.info.get("ONT_ADF")
        ont_adr = rec.info.get("ONT_ADR")
        if ont_adf and ont_adr:
            ref = ont_adf[0] + ont_adr[0]
            alt = ont_adf[1] + ont_adr[1]
            if ref + alt > 0:
                rec.info["ONT_VAF"] = float(alt / (ref + alt))
    except Exception:
        pass

    out.write(rec)

out.close()

print(f"Annotated VCF written to {args.out}")

