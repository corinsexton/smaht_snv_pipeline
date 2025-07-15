#!/usr/bin/env nextflow

/* 
 * Script to run VEP on split VCF files
SegDup="/n/data1/hms/dbmi/park/SOFTWARE/UCSC/GRCh38_UCSC_SegDup.formatted.bed"
 */

nextflow.enable.dsl=2

out_file = null
process filter_vep {
  cpus 1
  memory '2G'
  time '30m'

  label 'vep'

  publishDir "${params.results_dir}/vep_filtered",
    pattern: "${out_file}.gz*",
    mode: 'move'

  tag "$vcf.baseName"

  input:
    tuple(path(vcf), path(vcf_index))

  output:
    tuple( path("${out_file}.gz"), path("${out_file}.gz.tbi"))


  script:
  out_file = file(vcf).getSimpleName() + "_filtered.vcf"

    """
    filter_vep --input_file ${vcf} \
      --gz \
      --output_file ${out_file} \
      --force_overwrite --only_matched \
      -f "(MAX_AF < 0.001 or not MAX_AF) and (gnomad4.1_AF_grpmax_joint < 0.001 or not gnomad4.1_AF_grpmax_joint)"

    # Index the output
    bgzip ${out_file}
    tabix ${out_file}.gz
    """
}

