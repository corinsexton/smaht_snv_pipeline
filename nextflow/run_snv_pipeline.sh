#!/bin/bash


#SBATCH --job-name=nf_test
#SBATCH -A park
#SBATCH --partition short,park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 05:00:00
#SBATCH -o slurm-%x.%j.out


nextflow run main.nf -resume \
  --vep_config vep.ini \
  --input_csv truth_tests2.csv \
  --results_dir ./truthset_results

#nextflow run main.nf -resume \
#  --vep_config vep.ini \
#  --input_csv production_donors.csv \
#  --results_dir ./production_results
