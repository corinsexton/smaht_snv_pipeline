#!/bin/bash


#SBATCH --job-name=nf_benchmark
#SBATCH -A park
#SBATCH --partition park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -o slurm-%x.%j.out
##SBATCH --dependency=afterany:22694505


nextflow run main.nf  \
  --vep_config vep.ini \
  --longread_csv NEW_branch_truth_sets_lr_pooled.csv \
  --ont_csv NEW_branch_truth_sets_ont_pooled.csv \
  --shortread_csv NEW_branch_truth_sets_sr.csv \
  --input_metadata NEW_branch_truth_sets_metadata_PTA.csv \
  --input_vcfs NEW_branch_truth_sets_vcfs.csv \
  --results_dir ./NEW_branch_results_PTA_pileup
#
#
#nextflow run main.nf  -resume \
#  --vep_config vep.ini \
#  --longread_csv NEW_branch_truth_sets_lr_pooled.csv \
#  --ont_csv NEW_branch_truth_sets_ont_pooled.csv \
#  --shortread_csv NEW_branch_truth_sets_sr.csv \
#  --input_metadata NEW_branch_truth_sets_metadata_JINY.csv \
#  --input_vcfs NEW_branch_truth_sets_vcfs.csv \
#  --results_dir ./NEW_branch_results_JINY

#nextflow run main.nf  -resume \
#  --vep_config vep.ini \
#  --longread_csv Production_lr.csv \
#  --ont_csv Production_ont.csv \
#  --shortread_csv Production_sr.csv \
#  --input_metadata Production_metadata.csv \
#  --input_vcfs Production_vcfs.csv \
#  --results_dir ./NEW_branch_production_results


