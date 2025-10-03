process tier_variants {

    publishDir "${params.results_dir}/tiered_variants",
    pattern: "${id}.tiered.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/minipileup/parsed",
    pattern: "${id}.parsed.minipileup.tsv",
    mode:'copy'

    publishDir "${params.results_dir}/tiered_variants",
    pattern: "${id}.tier.metrics.tsv",
    mode:'copy'


    cache false

    cpus 1
    memory '1G'
    time '1h'

    tag "$id"

    input:
    tuple val(id), 
          path(vcf), path(tbi),
          path(truth_vcf), path(truth_vcf_tbi),
          path(minipileup_vcf)
    path(labels)

    output:
    tuple val(id), path("${id}.tiered.vcf.gz"), path("${id}.tiered.vcf.gz.tbi"),path(truth_vcf), path(truth_vcf_tbi), emit: vcf
    path("${id}.tier.metrics.tsv"), emit: metrics

    script:
    """
    labs=\$( cat ${labels} )

    # get the counts for pileups
    # matches ref and alt alleles
    parse_minipileup.py ${vcf} \
            ${minipileup_vcf} \
            ${id}.parsed.minipileup.tsv \
            --labels \${labs}

    # assigns tiers based on LR ≥ 1 (TIER1), LR = 0 & SR ≥ 2 (TIER2), else (.)
    split_lr_presence.py ${id}.parsed.minipileup.tsv ${id}.tiered.vcf \
        --original_vcf ${vcf}
    bgzip ${id}.tiered.vcf
    tabix ${id}.tiered.vcf.gz

    # split by tiers
    bcftools view -i 'FILTER="TIER1"' ${id}.tiered.vcf.gz -Oz -o ${id}.tier1.vcf.gz
    tabix ${id}.tier1.vcf.gz

    bcftools view -i 'FILTER="TIER2"' ${id}.tiered.vcf.gz -Oz -o ${id}.tier2.vcf.gz
    tabix ${id}.tier2.vcf.gz

    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    #  AFTER_VCF=${id}.tiered.vcf.gz
    TIER1_vcf=${id}.tier1.vcf.gz
    TIER2_vcf=${id}.tier2.vcf.gz

    num_before=\$(bcftools view -H "\${BEFORE_VCF}" | wc -l | awk '{print \$1}')
    # num_after=\$(bcftools view -H "\${AFTER_VCF}"  | wc -l | awk '{print \$1}')

    num_after_tier1=\$(bcftools view -H "\${TIER1_vcf}"  | wc -l | awk '{print \$1}')
    num_after_tier2=\$(bcftools view -H "\${TIER2_vcf}"  | wc -l | awk '{print \$1}')

    # Compute truth overlaps only if truth files are present
    if [[ -f "${truth_vcf}" ]]; then
      num_truth_before=\$(bcftools isec -n=2 -w1 -c both "\${BEFORE_VCF}" "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
      #num_truth_after=\$( bcftools isec -n=2 -w1 -c both "\${AFTER_VCF}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')

      num_truth_after_tier1=\$( bcftools isec -n=2 -w1 -c both "\${TIER1_vcf}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')

      num_truth_after_tier2=\$( bcftools isec -n=2 -w1 -c both "\${TIER2_vcf}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
    else
      num_truth_before=NA
      num_truth_after=NA
    fi

    {
      echo -e "id\tstep\tnum_before\tnum_truth_before\tnum_after\tnum_truth_after"
      # echo -e "${id}\tfilter_binom_fisher\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
      echo -e "${id}\tfiltered_tier1\t\${num_before}\t\${num_truth_before}\t\${num_after_tier1}\t\${num_truth_after_tier1}"
      echo -e "${id}\tfiltered_tier2\t\${num_before}\t\${num_truth_before}\t\${num_after_tier2}\t\${num_truth_after_tier2}"
    } > ${id}.tier.metrics.tsv


    """
}
