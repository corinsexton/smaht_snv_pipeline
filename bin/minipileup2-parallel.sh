#!/bin/bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -r reference.fasta [-o prefix] \\
          [--sr-cram CRAM ...] [--lr-cram CRAM ... --lr-tissue TISSUE --lr-type PB|ONT ...]

  -i             Input VCF (bgzipped) with .tbi index (required)
  -r             Reference FASTA with .fai index (required)
  -o             Output prefix (default: output)

  --sr-cram      Short-read CRAM with .crai index (repeatable)
  --lr-cram      Long-read CRAM (repeatable)
  --lr-tissue    Tissue ID matching each --lr-cram (repeatable)
  --lr-type      Sequencing type matching each --lr-cram (repeatable, PB or ONT)
EOF
  exit 1
}

INPUT_VCF=""
REFERENCE_FASTA=""
OUTPUT_PRFX="output"

SR_CRAMS=()
LR_CRAMS=()
LR_TISSUES=()
LR_TYPES=()
PB_CRAMS=()
PB_TISSUES=()
ONT_CRAMS=()
ONT_TISSUES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_VCF="$2"; shift 2;;
    -r) REFERENCE_FASTA="$2"; shift 2;;
    -o) OUTPUT_PRFX="$2"; shift 2;;
    --sr-cram)   SR_CRAMS+=("$2"); shift 2;;
    --lr-cram)   LR_CRAMS+=("$2"); shift 2;;
    --lr-tissue) LR_TISSUES+=("$2"); shift 2;;
    --lr-type)   LR_TYPES+=("$2"); shift 2;;
    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

[[ -n "$INPUT_VCF" ]]        || { echo "Error: -i required"; usage; }
[[ -n "$REFERENCE_FASTA" ]]  || { echo "Error: -r required"; usage; }

[[ -f "$INPUT_VCF" ]]             || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]]       || { echo "Error: ${INPUT_VCF}.tbi not found"; exit 1; }
[[ -f "$REFERENCE_FASTA" ]]       || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

(( ${#LR_CRAMS[@]} == ${#LR_TISSUES[@]} )) || { echo "Error: --lr-cram/--lr-tissue count mismatch"; exit 1; }
(( ${#LR_CRAMS[@]} == ${#LR_TYPES[@]} ))   || { echo "Error: --lr-cram/--lr-type count mismatch"; exit 1; }

for i in "${!LR_CRAMS[@]}"; do
  t_up="$(printf "%s" "${LR_TYPES[$i]}" | tr '[:lower:]' '[:upper:]')"
  case "$t_up" in
    PB)  PB_CRAMS+=("${LR_CRAMS[$i]}");  PB_TISSUES+=("${LR_TISSUES[$i]}");;
    ONT) ONT_CRAMS+=("${LR_CRAMS[$i]}"); ONT_TISSUES+=("${LR_TISSUES[$i]}");;
    *) echo "Error: invalid --lr-type '${LR_TYPES[$i]}' for '${LR_CRAMS[$i]}'. Expected PB or ONT."; exit 1;;
  esac
done

ALL_CRAMS=("${SR_CRAMS[@]}" "${PB_CRAMS[@]}" "${ONT_CRAMS[@]}")
(( ${#ALL_CRAMS[@]} > 0 )) || { echo "Error: no CRAM/BAM files provided"; exit 1; }

for b in "${ALL_CRAMS[@]}"; do
  [[ -f "$b" ]] || { echo "Error: file not found: $b"; exit 1; }
  if [[ "$b" == *.cram ]]; then
    [[ -f "${b}.crai" ]] || { echo "Error: missing .crai for $b"; exit 1; }
  elif [[ "$b" == *.bam ]]; then
    [[ -f "${b}.bai" || -f "${b}.csi" ]] || { echo "Error: missing BAM index for $b"; exit 1; }
  else
    echo "Error: unsupported file type (must be .bam or .cram): $b"; exit 1
  fi
done

command -v minipileup2 >/dev/null 2>&1 || { echo "Error: minipileup2 not in PATH"; exit 1; }
command -v bcftools    >/dev/null 2>&1 || { echo "Error: bcftools not in PATH"; exit 1; }
command -v bgzip       >/dev/null 2>&1 || { echo "Error: bgzip not in PATH"; exit 1; }
command -v tabix       >/dev/null 2>&1 || { echo "Error: tabix not in PATH"; exit 1; }
command -v python3     >/dev/null 2>&1 || { echo "Error: python3 not in PATH"; exit 1; }
python3 - <<'PY' || { echo "Error: granite not importable"; exit 1; }
from granite.lib import vcf_parser
PY

WORKDIR="$(mktemp -d minipileup2_work.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

echo "Running minipileup2..."
minipileup2 -t 2 -f "$REFERENCE_FASTA" -x "$INPUT_VCF" \
  -D -c -C -Q 20 -q 30 -s 0 \
  "${ALL_CRAMS[@]}" > "$WORKDIR/merged.vcf" \
  || { echo "Error: minipileup2 failed"; exit 1; }

[[ -s "$WORKDIR/merged.vcf" ]] || { echo "Error: minipileup2 produced empty output"; exit 1; }

bcftools sort -T "tmp_bcftools.XXXXXX" -O v -o "$WORKDIR/sorted.vcf" "$WORKDIR/merged.vcf" \
  || { echo "Error: bcftools sort failed"; exit 1; }

# Build sample map: stem -> type[-tissue]
: > "$WORKDIR/map.txt"
for cram in "${SR_CRAMS[@]}"; do
  printf "%s\tSR\t.\n" "$(basename "${cram%.*}")" >> "$WORKDIR/map.txt"
done
for i in "${!PB_CRAMS[@]}"; do
  printf "%s\tPB\t%s\n" "$(basename "${PB_CRAMS[$i]%.*}")" "${PB_TISSUES[$i]}" >> "$WORKDIR/map.txt"
done
for i in "${!ONT_CRAMS[@]}"; do
  printf "%s\tONT\t%s\n" "$(basename "${ONT_CRAMS[$i]%.*}")" "${ONT_TISSUES[$i]}" >> "$WORKDIR/map.txt"
done

py_script="
from granite.lib import vcf_parser
import os, re

sample_map = {}
with open('$WORKDIR/map.txt') as f:
    for line in f:
        stem, ftype, tissue = line.rstrip('\n').split('\t')
        sample_map[stem] = ftype if tissue == '.' else f'{ftype}-{tissue}'

vcf_obj = vcf_parser.Vcf('$WORKDIR/sorted.vcf')
cols = vcf_obj.header.columns.rstrip('\n').split('\t')
fixed = cols[:9]

for col in cols[9:]:
    stem = re.sub(r'\.(cram|bam)$', '', os.path.basename(col))
    ftype = sample_map[stem]
    fixed.append(f'{stem}-{ftype}')

vcf_obj.header.columns = '\t'.join(fixed) + '\n'

with open('${OUTPUT_PRFX}.vcf', 'w') as fo:
    vcf_obj.write_header(fo)
    for vnt in vcf_obj.parse_variants():
        vcf_obj.write_variant(fo, vnt)
"

echo "Renaming samples..."
python3 -c "$py_script" || { echo "Error: sample renaming failed"; exit 1; }

bgzip -f "${OUTPUT_PRFX}.vcf"   || { echo "Error: bgzip failed"; exit 1; }
tabix -f -p vcf "${OUTPUT_PRFX}.vcf.gz" || { echo "Error: tabix failed"; exit 1; }

echo "Done: ${OUTPUT_PRFX}.vcf.gz"
