#!/bin/bash


#SBATCH --job-name=SMHT004-3C_001_segdup
#SBATCH -A park
#SBATCH --partition park,short
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -o slurm-%x.%j.out
#SBATCH --dependency=afterany:30984163

#nextflow run main.nf  \
#  --vep_config vep.ini \
#  --longread_csv Production_lr.csv \
#  --ont_csv Production_ont.csv \
#  --shortread_csv Production_sr.csv \
#  --input_metadata Production_metadata.csv \
#  --input_vcfs Production_vcfs.csv \
#  --results_dir ./Production_results



#nextflow run main.nf -resume \
#  --vep_config vep.ini \
#  --longread_csv Production_lr.csv \
#  --ont_csv Production_ont.csv \
#  --shortread_csv Production_sr.csv \
#  --input_metadata SMHT005-3AK_metadata.csv \
#  --input_vcfs SMHT005-3AK_vcfs.csv \
#  --results_dir ./SMHT005-3AK_test_results


#nextflow run main.nf \
#  --vep_config vep.ini \
#  --longread_csv Production_lr.csv \
#  --ont_csv Production_ont.csv \
#  --shortread_csv Production_sr.csv \
#  --input_metadata SMHT005-3C_metadata.csv \
#  --input_vcfs SMHT005-3C_vcfs_norufus.csv \
#  --results_dir ./SMHT005-3C_test_results_noRUFUS/

#nextflow run main.nf -resume \
#  --vep_config vep.ini \
#  --longread_csv Production_lr.csv \
#  --ont_csv Production_ont.csv \
#  --shortread_csv Production_sr.csv \
#  --input_metadata Production_metadata.csv \
#  --input_vcfs SMHT004-3A_vcfs.csv \
#  --results_dir ./SMHT004-3A_test_results/

nextflow run main.nf \
  --vep_config vep.ini \
  --longread_csv Production_p5_lr.csv \
  --ont_csv Production_p5_ont.csv \
  --shortread_csv Production_p5_sr.csv \
  --input_metadata Production_p5_metadata.csv \
  --input_vcfs SMHT004-3C_vcfs.csv \
  --results_dir ./SMHT004-3C_test_results_segdup_1e-5
