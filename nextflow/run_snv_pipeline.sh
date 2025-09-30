#!/bin/bash


#SBATCH --job-name=nf_test
#SBATCH -A park
#SBATCH --partition short,park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 05:00:00
#SBATCH -o slurm-%x.%j.out


nextflow run main.nf  \
  --vep_config vep.ini \
  --longread_csv truth_sets_lr.csv \
  --ont_csv truth_sets_ont.csv \
  --shortread_csv truth_sets_sr.csv \
  --input_csv truth_sets_small.csv \
  --results_dir ./truthset_new_ref_files

#nextflow run main.nf -resume \
#  --vep_config vep.ini \
#  --longread_csv production_donors_lr.csv \
#  --ont_csv production_donors_ont.csv \
#  --shortread_csv production_donors_sr.csv \
#  --input_csv production_donors_small.csv \
#  --results_dir ./production_nosr_vaf_results
