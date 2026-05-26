#!/bin/bash


#SBATCH --job-name=nf_v2_ALL
#SBATCH -A park_contrib
#SBATCH --partition park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 7-00:00:00
#SBATCH -o logs/resume_NEW_v2-slurm-%x.%j.out
##SBATCH --dependency=afterany:38303592

nextflow run main.nf -resume \
  --vep_config vep.ini \
  --shortread_csv samplesheets/p25_sr.csv \
  --input_metadata samplesheets/p25_metadata.csv \
  --input_vcfs samplesheets/p25_vcfs.csv \
  --results_dir ./v2_p25_NOLR


  #--input_vcfs samplesheets/p25_vcfs.csv \
