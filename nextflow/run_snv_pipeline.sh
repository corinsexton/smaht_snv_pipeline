#!/bin/bash


#SBATCH --job-name=nf_PROD
#SBATCH -A park
#SBATCH --partition park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 36:00:00
#SBATCH -o slurm-%x.%j.out
#SBATCH --dependency=afterany:25861666

nextflow run main.nf  \
  --vep_config vep.ini \
  --longread_csv Production_lr.csv \
  --ont_csv Production_ont.csv \
  --shortread_csv Production_sr.csv \
  --input_metadata Production_metadata.csv \
  --input_vcfs Production_vcfs.csv \
  --results_dir ./Production_results
