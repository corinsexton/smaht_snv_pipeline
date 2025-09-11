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
    mode: 'copy'

  publishDir "${params.results_dir}/vep_filtered",
    pattern: "${id}.filter_vep.metrics.tsv",
    mode: 'copy'

  tag "$vcf.baseName"

  input:
    tuple(val(id), path(vcf), path(vcf_index), path(truth_vcf), path(truth_tbi))

  output:
    tuple val(id), path("${out_file}.gz"), path("${out_file}.gz.tbi"), path(truth_vcf), path(truth_tbi), emit: vcf
    path "${id}.filter_vep.metrics.tsv", emit: metrics
    


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


    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${out_file}.gz

    num_before=\$(bcftools view -H "\${BEFORE_VCF}" | wc -l | awk '{print \$1}')
    num_after=\$(bcftools view -H "\${AFTER_VCF}"  | wc -l | awk '{print \$1}')

    # Compute truth overlaps only if truth files are present
    if [[ -f "${truth_vcf}" ]]; then
      num_truth_before=\$(bcftools isec -n=2 -w1 -c both "\${BEFORE_VCF}" "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
      num_truth_after=\$( bcftools isec -n=2 -w1 -c both "\${AFTER_VCF}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
    else
      num_truth_before=NA
      num_truth_after=NA
    fi

    {
      echo -e "id\tstep\tnum_before\tnum_truth_before\tnum_after\tnum_truth_after"
      echo -e "${id}\tfilter_vep\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.filter_vep.metrics.tsv



    """
}

