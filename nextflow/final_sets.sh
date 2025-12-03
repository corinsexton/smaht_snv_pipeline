#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Usage info
# -------------------------------
usage() {
    echo "Usage: $0 -o output.vcf.gz input1.vcf.gz [input2.vcf.gz ...]"
    exit 1
}

# -------------------------------
# Parse arguments
# -------------------------------
if [[ $# -lt 3 ]]; then
    usage
fi

OUTFILE=""
VCFS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output)
            shift
            OUTFILE="${1:-}"
            [[ -z "$OUTFILE" ]] && usage
            shift
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            VCFS+=("$1")
            shift
            ;;
    esac
done

# -------------------------------
# Validate inputs
# -------------------------------
if [[ -z "$OUTFILE" ]]; then
    echo "Error: missing -o output file"
    usage
fi

if [[ ${#VCFS[@]} -eq 0 ]]; then
    echo "Error: no VCFs provided"
    usage
fi

for f in "${VCFS[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: file not found: $f"
        exit 1
    fi
done

# -------------------------------
# Run the Python script
# -------------------------------


./merge_callers.py -o tmp.vcf.gz "${VCFS[@]}"
tier1=$( zcat tmp.vcf.gz | grep 'TIER1' -c - )
tier2=$( zcat tmp.vcf.gz | grep 'TIER2' -c - )

./sandbox/remove_tier2_singletons.py -o "$OUTFILE" tmp.vcf.gz
tabix -f "$OUTFILE"
tier1_after=$( zcat "$OUTFILE" | grep 'TIER1' -c - )
tier2_after=$( zcat "$OUTFILE" | grep 'TIER2' -c - )

echo -e "$OUTFILE\tTier1 before\t${tier1}"
echo -e "$OUTFILE\tTier1 after\t${tier1_after}"
echo -e "$OUTFILE\tTier2 before\t${tier2}"
echo -e "$OUTFILE\tTier2 after\t${tier2_after}"
