#!/bin/bash


#SBATCH --job-name=nf_test
#SBATCH -A park
#SBATCH --partition short,park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 05:00:00
#SBATCH -o slurm-%x.%j.out


nextflow run main.nf -resume \
  --bin_size 10000 \
  --vep_config vep.ini \
  --input_csv 300x_files_w_bams.csv \
  --results_dir ./300x_results
