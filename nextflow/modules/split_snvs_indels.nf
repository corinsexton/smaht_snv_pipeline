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


  publishDir "${params.results_dir}/vep_filtered",
    pattern: "${id}.snv_only.metrics.tsv",
    mode: 'copy'

  tag "$vcf.baseName"

  input:
    tuple val(id), path(vcf), path(vcf_index), path(truth_vcf), path(truth_tbi)

  output:
    tuple val(id), path("snvs_${vcf.baseName}.gz"), path("snvs_${vcf.baseName}.gz.tbi"), emit: snvs
    tuple val(id), path("indels_${vcf.baseName}.gz"), path("indels_${vcf.baseName}.gz.tbi"), emit: indels

  script:
    """
    bcftools view -v snps ${vcf} --write-index=tbi -Ob -o snvs_${vcf.baseName}.gz
    bcftools view -v indels ${vcf}  --write-index=tbi -Ob -o indels_${vcf.baseName}.gz

    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=snvs_${vcf.baseName}.gz

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
      echo -e "${id}\tsnv_only\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.snv_only.metrics.tsv

    """
}

