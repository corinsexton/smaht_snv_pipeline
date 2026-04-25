#!/bin/bash
# Compare v2 multi-sample VCFs against v1 per-sample VCFs.
# Reports shared/unique variant counts and cross-tag breakdown per category.
#
# Usage:
#   ./compare_v1_v2.sh <v2_final_dir> <v1_filtered_calls_dir>
#
# v2_final_dir:          directory of *.final.vcf.gz files (one per tissue, e.g. SMHT001-3AD.final.vcf.gz)
# v1_filtered_calls_dir: directory organized as <DONOR>/*.vcf.gz, with tissue codes in filenames
#
# Example:
#   ./compare_v1_v2.sh ../results_test3/13_final ~/smaht/v1_filtered_calls

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <v2_final_dir> <v1_filtered_calls_dir>  || ./compare_v1_v2.sh ../results_test3/13_final ~/smaht/v1_filtered_calls" >&2
  exit 1
fi

V2_DIR="$1"
V1_DIR="$2"

TMPDIR_WORK=$(mktemp -d)
trap 'rm -rf "$TMPDIR_WORK"' EXIT

cross_tag_summary() {
  local vcf="$1"
  local coords="$2"
  local tags=("${@:3}")
  local n
  n=$(wc -l < "$coords")

  if [[ $n -eq 0 ]]; then
    for tag in "${tags[@]}"; do
      printf "    %s: 0/0\n" "$tag"
    done
    return
  fi

  local fmt=""
  for tag in "${tags[@]}"; do
    fmt+="%${tag}\t"
  done
  fmt="${fmt%\\t}"  # trim trailing \t

  bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t${fmt}\n" "$vcf" \
    | awk 'BEGIN{OFS="\t"} NR==FNR{key=$1"\t"$2"\t"$3"\t"$4; c[key]=1; next}
           {key=$1"\t"$2"\t"$3"\t"$4; if(key in c) print}' "$coords" - \
    | awk -v ntags="${#tags[@]}" -v tagnames="$(IFS=','; echo "${tags[*]}")" '
        BEGIN{ n=split(tagnames,tnames,",") }
        { for(i=1;i<=ntags;i++) { val=$(4+i); counts[i]+=(val=="1") } total++ }
        END{ for(i=1;i<=ntags;i++) printf "    %s: %d/%d\n", tnames[i], counts[i], total }'
}

filter_summary() {
  local vcf="$1"
  local coords="$2"
  local extra_field="${3:-}"
  local n
  n=$(wc -l < "$coords")

  if [[ $n -eq 0 ]]; then
    echo "    (none)"
    return
  fi

  local fmt="%CHROM\t%POS\t%REF\t%ALT\t%FILTER"
  [[ -n "$extra_field" ]] && fmt+="\t%${extra_field}"

  bcftools query -f "${fmt}\n" "$vcf" \
    | awk 'BEGIN{OFS="\t"} NR==FNR{key=$1"\t"$2"\t"$3"\t"$4; c[key]=1; next}
           {key=$1"\t"$2"\t"$3"\t"$4; if(key in c) {out=$5; if(NF>5) out=out" "$6; print out}}' \
        "$coords" - \
    | sort | uniq -c | sort -rn \
    | awk '{printf "    %s\t%s\n", $1, substr($0, index($0,$2))}'
}

for V2_VCF in "$V2_DIR"/*.final.vcf.gz; do
  BASENAME=$(basename "$V2_VCF" .final.vcf.gz)   # e.g. SMHT001-3AD
  DONOR=$(echo "$BASENAME" | cut -d- -f1)         # e.g. SMHT001
  TISSUE=$(echo "$BASENAME" | cut -d- -f2-)       # e.g. 3AD
  echo $TISSUE
  echo $DONOR
  # Find matching v1 VCF by tissue code in filename
  V1_VCF=$(find "$V1_DIR/$DONOR" -name "*.vcf.gz" ! -name "*.tbi" \
             -name "*-${TISSUE}-*" 2>/dev/null | head -1)
  echo $V1_VCF
  if [[ -z "$V1_VCF" ]]; then
    echo "=== $BASENAME === SKIPPED (no matching v1 VCF in $V1_DIR/$DONOR for tissue $TISSUE)"
    continue
  fi

  echo "=== $BASENAME ==="
  echo "  V2: $V2_VCF"
  echo "  V1: $V1_VCF"

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$V2_VCF" | sort > "$TMPDIR_WORK/v2.txt"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$V1_VCF" | sort > "$TMPDIR_WORK/v1.txt"

  comm -12 "$TMPDIR_WORK/v2.txt" "$TMPDIR_WORK/v1.txt" > "$TMPDIR_WORK/shared.txt"
  comm -23 "$TMPDIR_WORK/v2.txt" "$TMPDIR_WORK/v1.txt" > "$TMPDIR_WORK/unique_v2.txt"
  comm -13 "$TMPDIR_WORK/v2.txt" "$TMPDIR_WORK/v1.txt" > "$TMPDIR_WORK/unique_v1.txt"

  V2_TOTAL=$(wc -l < "$TMPDIR_WORK/v2.txt")
  V1_TOTAL=$(wc -l < "$TMPDIR_WORK/v1.txt")
  SHARED=$(wc -l < "$TMPDIR_WORK/shared.txt")
  UNIQUE_V2=$(wc -l < "$TMPDIR_WORK/unique_v2.txt")
  UNIQUE_V1=$(wc -l < "$TMPDIR_WORK/unique_v1.txt")

  printf "  V2 total: %d  |  V1 total: %d  |  Shared: %d  |  Unique V2: %d  |  Unique V1: %d\n" \
    "$V2_TOTAL" "$V1_TOTAL" "$SHARED" "$UNIQUE_V2" "$UNIQUE_V1"

  echo ""
  echo "  [SHARED] V2 cross tags:"
  cross_tag_summary "$V2_VCF" "$TMPDIR_WORK/shared.txt"    CrossTech CrossCore CrossTissue
  echo "  [SHARED] V2 filter/tier:"
  filter_summary    "$V2_VCF" "$TMPDIR_WORK/shared.txt"    TIER

  echo "  [UNIQUE V2] V2 cross tags:"
  cross_tag_summary "$V2_VCF" "$TMPDIR_WORK/unique_v2.txt" CrossTech CrossCore CrossTissue
  echo "  [UNIQUE V2] V2 filter/tier:"
  filter_summary    "$V2_VCF" "$TMPDIR_WORK/unique_v2.txt" TIER

  echo "  [SHARED] V1 cross tags:"
  cross_tag_summary "$V1_VCF" "$TMPDIR_WORK/shared.txt"    CrossTech CrossCaller CrossTissue
  echo "  [SHARED] V1 filter:"
  filter_summary    "$V1_VCF" "$TMPDIR_WORK/shared.txt"

  echo "  [UNIQUE V1] V1 cross tags:"
  cross_tag_summary "$V1_VCF" "$TMPDIR_WORK/unique_v1.txt" CrossTech CrossCaller CrossTissue
  echo "  [UNIQUE V1] V1 filter:"
  filter_summary    "$V1_VCF" "$TMPDIR_WORK/unique_v1.txt"

  echo ""
done
