#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------
# Step 4 — Map filtered germline VCF to TIER1 sites
#   Input: filtered germline VCF (AF > 0.001)
#   Output: final mapping table linking each TIER1 variant
#            to the *closest* germline variant (within ±5 kb)
# ---------------------------------------------------------------

usage() {
  echo "Usage: phasing_step4.sh -s sample -v pass.vcf.gz -w tier1_windows.bed"
}

while getopts ":s:v:w:" opt; do
  case $opt in
	s) SAMPLE=$OPTARG ;;
    v) PASS_VCF=$OPTARG ;;
    w) WIN_BED=$OPTARG ;;
    *) usage; exit 1 ;;
  esac
done


[[ -z ${PASS_VCF-} || -z ${WIN_BED-} ]] && { usage; exit 1; }

PASS_BED="${SAMPLE}.germline.pass.bed"
SOMATIC_BED="${SAMPLE}.somatic_sites.bed"
CLOSEST="${SAMPLE}.closest.tsv"
MAP_TSV="${SAMPLE}.germline_map.tsv"
FAILED_TSV="${SAMPLE}.no_close_germline.tsv"

echo "[Step4] Converting filtered germline VCF to BED..."
bcftools query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' "$PASS_VCF" > "$PASS_BED"

echo "[Step4] Generating 1bp BED for somatic sites..."
awk 'BEGIN{OFS="\t"}{print $1,$5-1,$5,$5,$6,$7}' "$WIN_BED" > "$SOMATIC_BED"

bedtools sort -i "$PASS_BED" > z; mv z "$PASS_BED"
bedtools sort -i "$SOMATIC_BED" > z; mv z "$SOMATIC_BED"

echo "[Step4] Finding closest germline variant to each somatic site..."
bedtools closest -a "$SOMATIC_BED" -b "$PASS_BED" -d | awk '$12 != -1 && $12 <= 5000' > "$CLOSEST"

# Handle somatic variants with no germline within 5kb
echo "[Step4] Identifying somatic sites with no germline within 5kb..."
bedtools closest -a "$SOMATIC_BED" -b "$PASS_BED" -d | awk '$12 == -1 || $12 > 5000' > "$FAILED_TSV"

# Process passing (within 5kb) and failing (>5kb) separately, then combine
awk '
BEGIN {OFS="\t"}
{
  var_chrom=$1; var_pos=$4; var_ref=$5; var_alt=$6;
  germ_pos=$8+1; germ_ref=$10; germ_alt=$11;
  if (germ_ref == "" || germ_ref == ".") {
    germ_pos="NA"; germ_ref="NA"; germ_alt="NA";
  }
  print var_chrom, var_pos, var_ref, var_alt, germ_pos, germ_ref, germ_alt;
}' "$CLOSEST" | sort -k1,1 -k2,2n > tmp.pass.tsv

awk '
BEGIN {OFS="\t"}
{
  var_chrom=$1; var_pos=$4; var_ref=$5; var_alt=$6;
  print var_chrom, var_pos, var_ref, var_alt, "NA", "NA", "NA";
}' "$FAILED_TSV" | sort -k1,1 -k2,2n > tmp.fail.tsv

{
  echo -e "chrom\tVar_pos\tVar_ref\tVar_alt\tGerm_pos\tGerm_ref\tGerm_alt"
  cat tmp.pass.tsv tmp.fail.tsv
} > "$MAP_TSV"

rm -f tmp.pass.tsv tmp.fail.tsv

echo "[Step4 complete] Final mapping written to $MAP_TSV (within ±5kb only)"

