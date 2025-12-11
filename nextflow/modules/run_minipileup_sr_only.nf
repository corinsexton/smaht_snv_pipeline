process run_minipileup_sr_only {

    publishDir "${params.results_dir}/minipileup_sr",
    pattern: "${id}.minipileup_sr.vcf.gz",
    mode:'copy'

    publishDir "${params.results_dir}/final",
    pattern: "${id}.final.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/final",
    pattern: "${id}.final.*.tsv",
    mode:'copy'

    cache 'lenient'


    cpus 4 
    memory '4G'
    time '6h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), 
        path(truth_vcf), path(truth_vcf_tbi),
        path(sr_bams), path(sr_bais), val(sr_ids)
    tuple path(ref), path(ref_index), path(ref_dict)
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)


    output:
    tuple val(id),
          path("${id}.final.vcf.gz"), path("${id}.final.vcf.gz.tbi"),
          path(truth_vcf), path(truth_vcf_tbi), emit: vcf
    path("${id}.minipileup_sr.vcf.gz")
    path("${id}.final.metrics.tsv")
    path("${id}.final.regions.tsv")

    script:
    """
    current_tissue=\$(echo "$id" | cut -d'-' -f2)

    # Build: --sr-cram <bam1> --sr-cram <bam2> ...
    sr_crams=""
    for f in ${sr_bams}; do
        sr_crams+=" --sr-cram \${f}"
    done

    # Build: --sr-tissue <tissue id1> --sr-tissue <tissue id2> ...
    sr_tissue=""
    for f in ${sr_ids}; do
        sr_tissue+=" --sr-tissue \${f}"
    done

    minipileup-parallel_sr_only.sh -i ${vcf} \
        -r ${ref} \
        -t 2 \
        -o ${id}.minipileup_sr \
        \${sr_crams} \
        \${sr_tissue}

    parse_minipileup_sr_only.py \
        --tissue \${current_tissue} \
        --orig_vcf ${vcf} \
        --mp_vcf ${id}.minipileup_sr.vcf.gz \
        --out ${id}.CrossTissue.vcf.gz

    tabix ${id}.CrossTissue.vcf.gz

    set_filter_vcf.py \
        -i ${id}.CrossTissue.vcf.gz \
        -o ${id}.final.vcf.gz

    tabix ${id}.final.vcf.gz

    # split by tiers
    bcftools view -i 'FILTER="HighConf"' ${id}.final.vcf.gz -Oz -o ${id}.tier1.vcf.gz
    tabix ${id}.tier1.vcf.gz

    bcftools view -i 'FILTER="ModerateConf"' ${id}.final.vcf.gz -Oz -o ${id}.tier2.vcf.gz
    tabix ${id}.tier2.vcf.gz

    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.final.vcf.gz
    TIER1_vcf=${id}.tier1.vcf.gz
    TIER2_vcf=${id}.tier2.vcf.gz

    # check regions
    num_easy_before=\$( bedtools intersect -u -b $easy_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_easy_after=\$( bedtools intersect -u -b $easy_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_diff_before=\$( bedtools intersect -u -b $diff_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_diff_after=\$( bedtools intersect -u -b $diff_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_ext_before=\$( bedtools intersect -u -b $ext_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_ext_after=\$( bedtools intersect -u -b $ext_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')


    {
      echo -e "id\tstep\tregion_type\tnum_before\tnum_after"
      echo -e "${id}\tfinal\teasy\t\${num_easy_before}\t\${num_easy_after}"
      echo -e "${id}\tfinal\tdiff\t\${num_diff_before}\t\${num_diff_after}"
      echo -e "${id}\tfinal\text\t\${num_ext_before}\t\${num_ext_after}"
    } > ${id}.final.regions.tsv


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
      # echo -e "${id}\tfinal\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
      echo -e "${id}\tfinalhighconf\t\${num_before}\t\${num_truth_before}\t\${num_after_tier1}\t\${num_truth_after_tier1}"
      echo -e "${id}\tfinal_modconf\t\${num_before}\t\${num_truth_before}\t\${num_after_tier2}\t\${num_truth_after_tier2}"
    } > ${id}.final.metrics.tsv


    """
}
