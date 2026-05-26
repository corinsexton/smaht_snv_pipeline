#!/bin/bash

#SBATCH --job-name=nf_v2_test_nolr
#SBATCH -A park_contrib
#SBATCH --partition park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/test_nolr-slurm-%x.%j.out

# SR-only test run: SMHT001-3A and SMHT001-3I only (no LR CSVs passed)
nextflow run main.nf -resume \
  --vep_config vep.ini \
  --shortread_csv samplesheets/p25_sr.csv \
  --input_metadata samplesheets/p25_metadata.csv \
  --input_vcfs samplesheets/p25_vcfs_test.csv \
  --results_dir ./v2_test_NOLR
