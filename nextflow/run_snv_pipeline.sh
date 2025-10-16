#!/bin/bash


#SBATCH --job-name=nf_truth_sets
#SBATCH -A park
#SBATCH --partition short,park
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -o slurm-%x.%j.out


#nextflow run main.nf \
#  --vep_config vep.ini \
#  --longread_csv truth_sets_lr_pooled.csv \
#  --ont_csv truth_sets_ont_pooled.csv \
#  --shortread_csv truth_sets_sr.csv \
#  --input_csv truth_sets_small.csv \
#  --results_dir ./pooled_truthset_results
#
#nextflow run main.nf -resume \
#  --vep_config vep.ini \
#  --longread_csv truth_sets_lr_pooled.csv \
#  --ont_csv truth_sets_ont_pooled.csv \
#  --shortread_csv truth_sets_sr.csv \
#  --input_csv truth_sets_small_PTA.csv \
#  --results_dir ./pooled_truthset_PTA_results


  #--input_csv truth_sets_small.csv \
  #--results_dir ./truthset_results 

  #--input_csv truth_sets_small.csv \
  #--results_dir ./truthset_results 


  #--input_csv truth_sets_small_PTA_norufus.csv \
  #--results_dir ./truthset_PTA_results



nextflow run main.nf -resume \
  --vep_config vep.ini \
  --longread_csv production_donors_lr.csv \
  --ont_csv production_donors_ont.csv \
  --shortread_csv production_donors_sr.csv \
  --input_csv production_donors_small.csv \
  --results_dir ./production_results

nextflow run main.nf -resume \
  --vep_config vep.ini \
  --longread_csv prod_donors_150x_lr.csv \
  --ont_csv prod_donors_150x_ont.csv \
  --shortread_csv prod_donors_150x_sr.csv \
  --input_csv prod_donors_150x.csv \
  --results_dir ./prod_results_150x



#nextflow run main.nf -resume \
#  --vep_config vep.ini \
#  --longread_csv colo_lr.csv \
#  --shortread_csv colo_sr.csv \
#  --ont_csv colo_ont.csv \
#  --input_csv colo_small.csv \
#  --results_dir ./colo_results
