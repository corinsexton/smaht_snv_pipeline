#!/bin/bash


#SBATCH --job-name=nf_p25_big_vcf
#SBATCH -A park_contrib
#SBATCH --partition park
#SBATCH --mem 8G
#SBATCH -c 1
#SBATCH -t 30-00:00:00
#SBATCH -o slurm-%x.%j.out
##SBATCH --dependency=afterany:25861666

#nextflow run main.nf -resume \
#  --vep_config vep.ini \
#  --longread_csv samplesheets/p25/p25_lr.csv \
#  --ont_csv samplesheets/p25/p25_ont.csv \
#  --shortread_csv samplesheets/p25/p25_sr.csv \
#  --input_metadata samplesheets/p25/p25_metadata.csv \
#  --input_vcfs samplesheets/p25/p25_vcfs.csv \
#  --results_dir results_p25_big_vcf/


  #--input_vcfs p25_samplesheets_merged/p25_vcfs.csv \
  #--results_dir ./all_filter_results_p25

nextflow run main.nf -resume \
  --vep_config vep.ini \
  --longread_csv samplesheets/p25/p25_lr.csv \
  --ont_csv samplesheets/p25/p25_ont.csv \
  --shortread_csv samplesheets/p25/p25_sr.csv \
  --input_metadata samplesheets/p25/p25_metadata.csv \
  --input_vcfs samplesheets/p25/p25_vcfs_subset.csv \
  --results_dir results_p25_big_vcf/


