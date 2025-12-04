#!/bin/bash


#SBATCH --job-name=nf_production
#SBATCH -A park
#SBATCH --partition park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -o slurm-%x.%j.out
#SBATCH --dependency=afterany:22694505


nextflow run main.nf  -resume \
  --vep_config vep.ini \
  --longread_csv NEW_branch_truth_sets_lr_pooled.csv \
  --ont_csv NEW_branch_truth_sets_ont_pooled.csv \
  --shortread_csv NEW_branch_truth_sets_sr.csv \
  --input_metadata NEW_branch_truth_sets_metadata.csv \
  --input_vcfs NEW_branch_truth_sets_vcfs.csv \
  --results_dir ./NEW_branch_results

