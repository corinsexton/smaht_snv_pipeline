#!/bin/bash


#SBATCH --job-name=nf_skin
#SBATCH -A park_contrib
#SBATCH --partition priopark
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 96:00:00
#SBATCH -o logs/slurm-%x.%j.out

nextflow run main.nf  -resume \
  --vep_config vep.ini \
  --longread_csv samplesheets_fixed/p25_lr.csv \
  --ont_csv samplesheets_fixed/p25_ont.csv \
  --shortread_csv samplesheets_fixed/p25_sr.csv \
  --input_metadata samplesheets/p25_metadata.csv \
  --input_vcfs p25_vcfs_ss_3AF.csv \
  --results_dir ./results_skin_mp2
