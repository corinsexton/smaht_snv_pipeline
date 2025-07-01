#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Usage: ./script.sh input.vcf output_directory
input_vcf=$1
DIR_vcf=$2

# Extract sample name (portion before first dot)
sample_name=$(basename "$input_vcf")
sample_name=${sample_name%%.*}

# Your BED directories
DIR_bed1=/n/data1/hms/dbmi/park/clara_kim/smaht/truth_set_generation/regions/051925_hengs/pm151_genomes
DIR_bed2=/n/data1/hms/dbmi/park/clara_kim/smaht/truth_set_generation/regions/061324_regions_to_distribute
BED_DIRS=("$DIR_bed1" "$DIR_bed2")

# Make sure output dir exists
mkdir -p "$DIR_vcf"

# Prepare bcftools path
bcftools=/n/data1/hms/dbmi/park/dominika/smaht_env/smaht_domini/bin/bcftools

# Summary file
summary_file="$DIR_vcf/${sample_name}_region_counts.tsv"
echo -e "Sample	Suffix	Count" > "$summary_file"

# Count total variants (optional if you want to include a total line)
total_variants=$($bcftools view -H "$input_vcf" | wc -l)
echo -e "${sample_name}	TOTAL	${total_variants}" >> "$summary_file"

# Process each BED
for bed_dir in "${BED_DIRS[@]}"; do
  beds=("$bed_dir"/*.bed)
  if [ ${#beds[@]} -eq 0 ]; then
    continue
  fi

  for bed in "${beds[@]}"; do
    region=$(basename "$bed" .bed)
    out_vcf="$DIR_vcf/${sample_name}.${region}.vcf"

    # Generate per-region VCF and count
    bedtools intersect -header -a "$input_vcf" -b "$bed" > "$out_vcf"
    count=$($bcftools view -H "$out_vcf" | wc -l)

    # Append to summary file
    echo -e "${sample_name}	${region}	${count}" >> "$summary_file"
  done
done

# Final message
echo "Done. Summary written to $summary_file"
