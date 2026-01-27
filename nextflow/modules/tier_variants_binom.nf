process tier_variants_binom {

    publishDir "${params.results_dir}/tiered_variants",
    pattern: "${id}.tiered.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/tiered_variants",
    pattern: "${id}.tier.metrics.tsv",
    mode:'copy'

    cpus 1
    memory '2G'
    time '1h'

    tag "$id"

    input:
    tuple val(id), 
          path(vcf), path(tbi),
          path(truth_vcf), path(truth_vcf_tbi),
          path(minipileup_vcf), path(minipileup_vcf_tbi)
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)

    output:
    tuple val(id), path("${id}.tiered.vcf.gz"), path("${id}.tiered.vcf.gz.tbi"),path(truth_vcf), path(truth_vcf_tbi), emit: vcf
    path("${id}.tier.metrics.tsv"), emit: metrics

    script:
    """


    tier_filter_variants_SR_PB_ONT.py \
        -i ${vcf} \
        -m ${minipileup_vcf} \
        --current_tissue ${id} \
        -o "${id}.tiered.vcf.gz" 


    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.tiered.vcf.gz

    # check regions
    num_easy_before=\$( bedtools intersect -u -b $easy_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_easy_after=\$( bedtools intersect -u -b $easy_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_diff_before=\$( bedtools intersect -u -b $diff_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_diff_after=\$( bedtools intersect -u -b $diff_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_ext_before=\$( bedtools intersect -u -b $ext_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_ext_after=\$( bedtools intersect -u -b $ext_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')


    {
      echo -e "id\tstep\tregion_type\tnum_before\tnum_after"
      echo -e "${id}\ttier\teasy\t\${num_easy_before}\t\${num_easy_after}"
      echo -e "${id}\ttier\tdiff\t\${num_diff_before}\t\${num_diff_after}"
      echo -e "${id}\ttier\text\t\${num_ext_before}\t\${num_ext_after}"
    } > ${id}.tier.regions.tsv


    num_before=\$(bcftools view -H "\${BEFORE_VCF}" | wc -l | awk '{print \$1}')
    num_after=\$(bcftools view -H "\${AFTER_VCF}"  | wc -l | awk '{print \$1}')

    # Compute truth overlaps only if truth files are present
    if [[ -f "${truth_vcf}" ]]; then
      num_truth_before=\$(bcftools isec -n=2 -w1 -c both "\${BEFORE_VCF}" "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
      num_truth_after=\$( bcftools isec -n=2 -w1 -c both "\${AFTER_VCF}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')

      #num_truth_after_tier1=\$( bcftools isec -n=2 -w1 -c both "\${TIER1_vcf}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')

      #num_truth_after_tier2=\$( bcftools isec -n=2 -w1 -c both "\${TIER2_vcf}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
    else
      num_truth_before=NA
      num_truth_after=NA
    fi

    {
      echo -e "id\tstep\tnum_before\tnum_truth_before\tnum_after\tnum_truth_after"
      echo -e "${id}\tfilter_binom_fisher\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.tier.metrics.tsv


    """
}
