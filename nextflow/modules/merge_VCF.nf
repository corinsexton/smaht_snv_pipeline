#!/usr/bin/env nextflow

/* 
 * Script to merge VCF files into a single file
 */

nextflow.enable.dsl=2

// defaults
merged_vcf = null
if ( params.output_prefix ){
  merged_vcf = params.output_prefix + "_VEP.vcf.gz"
}

process mergeVCF {
  /*
  Merge VCF files into a single file
  */
    
  publishDir "${params.results_dir}/vep_annotated",
    pattern: "${merged_vcf}*",
    mode: 'move'


  cpus params.cpus
  label 'bcftools'
  cache 'lenient'
   
  input:
  tuple val(id), val(original_file), path(vcf_files), path(index_files)
  
  output:
  path("${merged_vcf}*")

  script:

  index_type = "tbi"
  
  merged_vcf = merged_vcf ?: file(original_file).getSimpleName() + "_VEP.vcf.gz"
  index_flag = index_type == "tbi" ? "-t" : "-c"
  
  """
  sorted_vcfs=\$(echo ${vcf_files} | xargs -n1 | sort | xargs)
  bcftools concat --no-version --naive \${sorted_vcfs} -Oz -o ${merged_vcf}
  bcftools index ${index_flag} ${merged_vcf}
  """
}
