#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process split_snvs_indels {
  cpus 1
  memory '2G'
  time '30m'

  publishDir "${params.results_dir}/vep_filtered/snvs",
    pattern: "snvs_*.vcf.gz*",
    mode: 'copy'

  publishDir "${params.results_dir}/vep_filtered/indels",
    pattern: "indels_*.vcf.gz*",
    mode: 'copy'

  tag "$vcf.baseName"

  input:
    tuple val(id), path(vcf), path(vcf_index)

  output:
    tuple val(id), path("snvs_${vcf.baseName}.gz"), path("snvs_${vcf.baseName}.gz.tbi"), emit: snvs
    tuple val(id), path("indels_${vcf.baseName}.gz"), path("indels_${vcf.baseName}.gz.tbi"), emit: indels

  script:
    """
    bcftools view -v snps ${vcf} --write-index=tbi -Ob -o snvs_${vcf.baseName}.gz
    bcftools view -v indels ${vcf}  --write-index=tbi -Ob -o indels_${vcf.baseName}.gz
    """
}

