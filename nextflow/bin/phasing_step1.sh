#!/usr/bin/env bash
set -euo pipefail

usage() {
echo "Usage: phasing_step1_prep.sh -r ref -s sample -v somatic.vcf.gz -b bam_list [-t threads]"
}

THREADS=4
while getopts ":r:s:v:b:t:" opt; do
case $opt in
r) REF=$OPTARG;; s) SAMPLE=$OPTARG;; v) SOMATIC=$OPTARG;; b) BAM_LIST=$OPTARG;; t) THREADS=$OPTARG;; *) usage; exit 1;;
esac
done

LOG="step1.log"

TIER1_VCF="tier1.norm.vcf.gz"
WIN_BED="$SAMPLE.tier1_windows.bed"
WIN_BED_MERGED="tier1_windows.merged.bed"

bcftools view -i 'FILTER=="TIER1"' -Oz -o "$TIER1_VCF" "$SOMATIC"
tabix -f "$TIER1_VCF"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$TIER1_VCF" | \
awk 'BEGIN{OFS="\t"}{s=$2-5000;if(s<0)s=0;e=$2+5000;print $1,s,e,$1":"$2":"$3":"$4,$2,$3,$4}' > "$WIN_BED"

sort -k1,1 -k2,2n "$WIN_BED" | bedtools merge -i - -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "$WIN_BED_MERGED"

HC_VCF="$SAMPLE.hc.vcf.gz"

/n/data1/hms/dbmi/park/SOFTWARE/GATK/gatk-4.6.1.0/gatk --java-options "-Xmx12g" HaplotypeCaller \
	-R "$REF" -I ${BAM_LIST} -L "$WIN_BED_MERGED" --native-pair-hmm-threads "$THREADS" -O "$HC_VCF"

echo "[STEP1] Filtering to heterozygous SNVs..."
bcftools norm -m -both -f "$REF" "$HC_VCF" -Ou \
  | bcftools view -i 'TYPE="snp" && GT="het"' -Oz -o "$SAMPLE.hc.norm.vcf.gz"
tabix -f "$SAMPLE.hc.norm.vcf.gz"


#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$SAMPLE.hc.norm.vcf.gz" > "$SAMPLE.hc.norm.simple.tsv"

echo "[Step1 complete]"
