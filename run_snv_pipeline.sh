#!/bin/bash


#SBATCH --job-name=nf_big_vcf_test
#SBATCH -A park_contrib
#SBATCH --partition park
#SBATCH --mem 8G
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o slurm-%x.%j.out
##SBATCH --dependency=afterany:25861666

nextflow run main.nf -resume \
  --vep_config vep.ini \
  --longread_csv p25_lr_FIX.csv \
  --ont_csv p25_samplesheets_merged/p25_ont.csv \
  --shortread_csv p25_samplesheets_merged/p25_sr.csv \
  --input_metadata p25_samplesheets_merged/p25_metadata.csv \
  --input_vcfs p25_vcfs_ss.csv \
  --results_dir results_test


  #--input_vcfs p25_samplesheets_merged/p25_vcfs.csv \
  #--results_dir ./all_filter_results_p25
