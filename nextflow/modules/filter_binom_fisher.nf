process filter_binom_fisher {

    publishDir "${params.results_dir}/binom_fisher",
    pattern: "${id}.tiered.binom_fisher.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/binom_fisher",
    pattern: "${id}.filter_binom_fisher.metrics.tsv",
    mode:'copy'

    publishDir "${params.results_dir}/binom_fisher/failed_variants",
    pattern: "${id}.tiered.binom_fisher.vcf_failed_*.vcf",
    mode:'copy'

    cpus 1
    memory '1G'
    time '1h'

    tag "$id"

    input:
    tuple val(id), path(tiered_vcf), path(tiered_vcf_index), path(truth_vcf), path(truth_tbi)

    output:
    tuple val(id), path("${id}.tiered.binom_fisher.vcf.gz"),path("${id}.tiered.binom_fisher.vcf.gz.tbi"), path(truth_vcf), path(truth_tbi)
    path("${id}.filter_binom_fisher.metrics.tsv"), emit: metrics
    path("${id}.tiered.binom_fisher.vcf_failed_both.vcf")
    path("${id}.tiered.binom_fisher.vcf_failed_ont.vcf")
    path("${id}.tiered.binom_fisher.vcf_failed_pb.vcf")

    script:
    """
    filter_binom_fisher.py ${tiered_vcf} ${id}.tiered.binom_fisher.vcf
    bgzip ${id}.tiered.binom_fisher.vcf
    tabix ${id}.tiered.binom_fisher.vcf.gz


    # split by tiers
    bcftools view -i 'FILTER="TIER1"' ${id}.tiered.binom_fisher.vcf.gz -Oz -o ${id}.tier1.vcf.gz
    tabix ${id}.tier1.vcf.gz

    bcftools view -i 'FILTER="TIER2"' ${id}.tiered.binom_fisher.vcf.gz -Oz -o ${id}.tier2.vcf.gz
    tabix ${id}.tier2.vcf.gz

    # --- metrics (standard schema) ---
    BEFORE_VCF=${tiered_vcf}
    AFTER_VCF=${id}.tiered.binom_fisher.vcf.gz
    TIER1_vcf=${id}.tier1.vcf.gz
    TIER2_vcf=${id}.tier2.vcf.gz

    num_before=\$(bcftools view -H "\${BEFORE_VCF}" | wc -l | awk '{print \$1}')
    num_after=\$(bcftools view -H "\${AFTER_VCF}"  | wc -l | awk '{print \$1}')

    num_after_tier1=\$(bcftools view -H "\${TIER1_vcf}"  | wc -l | awk '{print \$1}')
    num_after_tier2=\$(bcftools view -H "\${TIER2_vcf}"  | wc -l | awk '{print \$1}')

    # Compute truth overlaps only if truth files are present
    if [[ -f "${truth_vcf}" ]]; then
      num_truth_before=\$(bcftools isec -n=2 -w1 -c both "\${BEFORE_VCF}" "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
      num_truth_after=\$( bcftools isec -n=2 -w1 -c both "\${AFTER_VCF}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')

      num_truth_after_tier1=\$( bcftools isec -n=2 -w1 -c both "\${TIER1_vcf}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')

      num_truth_after_tier2=\$( bcftools isec -n=2 -w1 -c both "\${TIER2_vcf}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
    else
      num_truth_before=NA
      num_truth_after=NA
    fi

    {
      echo -e "id\tstep\tnum_before\tnum_truth_before\tnum_after\tnum_truth_after"
      echo -e "${id}\tfilter_binom_fisher\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
      echo -e "${id}\tfilter_binom_fisher_tier1\t\${num_before}\t\${num_truth_before}\t\${num_after_tier1}\t\${num_truth_after_tier1}"
      echo -e "${id}\tfilter_binom_fisher_tier2\t\${num_before}\t\${num_truth_before}\t\${num_after_tier2}\t\${num_truth_after_tier2}"
    } > ${id}.filter_binom_fisher.metrics.tsv

    """
}
