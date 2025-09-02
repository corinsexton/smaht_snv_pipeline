#!/usr/bin/env nextflow

/* 
 * Script to run VEP on split VCF files
SegDup="/n/data1/hms/dbmi/park/SOFTWARE/UCSC/GRCh38_UCSC_SegDup.formatted.bed"
 */

nextflow.enable.dsl=2

out_file = null
process split_snvs_indels {
  cpus 1
  memory '2G'
  time '30m'

  publishDir "${params.results_dir}/vep_filtered/snvs",
    pattern: "snvs_*${out_file}.vcf.gz*",
    mode: 'move'

  publishDir "${params.results_dir}/vep_filtered/indels",
    pattern: "indels_*${out_file}.vcf.gz*",
    mode: 'move'

  tag "$vcf.baseName"

  input:
    tuple(path(vcf), path(vcf_index))

  output:
    tuple( path("snvs_${out_file}.vcf.gz"), path("snvs_${out_file}.vcf.gz.tbi"),
            path("indels_${out_file}.vcf.gz"), path("indels_${out_file}.vcf.gz.tbi"))


  script:
  out_file = file(vcf).getSimpleName()

    """
    bcftools view -v snps ${out_file}.vcf.gz --write-index=tbi -Ob -o snvs_${out_file}.vcf.gz
    bcftools view -v indels ${out_file}.vcf.gz --write-index=tbi -Ob -o indels_${out_file}.vcf.gz
    """
}

