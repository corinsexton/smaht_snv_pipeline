#!/bin/bash


#SBATCH --job-name=nf_big_vcf
#SBATCH -A park
#SBATCH --partition park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 36:00:00
#SBATCH -o slurm-%x.%j.out
##SBATCH --dependency=afterany:25861666

nextflow run main.nf -resume  \
  --vep_config vep.ini \
  --longread_csv p25_samplesheets_merged/p25_lr.csv \
  --ont_csv p25_samplesheets_merged/p25_ont.csv \
  --shortread_csv p25_samplesheets_merged/p25_sr.csv \
  --input_metadata p25_samplesheets_merged/p25_metadata.csv \
  --input_vcfs p25_vcfs_ss.csv \
  --results_dir ./results_p25_ss
