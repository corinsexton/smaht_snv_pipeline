process filter_germline_variants {


    publishDir "${params.results_dir}/germline_filtered",
    pattern: "${id}.nogermline.vcf.gz*",
    mode:'copy'


    publishDir "${params.results_dir}/germline_filtered",
    pattern: "${id}.filter_germline.*.tsv",
    mode:'copy'


    cpus 1
    memory '4G'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi),
        path(germline_vcf), path(germline_tbi)
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)
    

    output:
    tuple val(id),
          path("${id}.nogermline.vcf.gz"),
          path("${id}.nogermline.vcf.gz.tbi"),
          path(truth_vcf), path(truth_tbi), emit: vcf
        path("${id}.filter_germline.metrics.tsv"), emit: metrics
        path("${id}.filter_germline.regions.tsv"), emit: regions

    script:
    """

    bcftools isec -n~01 -p isec_tmp  \
        ${germline_vcf} \
        ${vcf}

    mv isec_tmp/0001.vcf ${id}.nogermline.vcf
    bgzip ${id}.nogermline.vcf
    tabix ${id}.nogermline.vcf.gz


    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.nogermline.vcf.gz

    # check regions
    num_easy_before=\$( bedtools intersect -u -b $easy_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_easy_after=\$( bedtools intersect -u -b $easy_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_diff_before=\$( bedtools intersect -u -b $diff_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_diff_after=\$( bedtools intersect -u -b $diff_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_ext_before=\$( bedtools intersect -u -b $ext_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_ext_after=\$( bedtools intersect -u -b $ext_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    {
      echo -e "id\tstep\tregion_type\tnum_before\tnum_after"
      echo -e "${id}\tfilter_germline\teasy\t\${num_easy_before}\t\${num_easy_after}"
      echo -e "${id}\tfilter_germline\tdiff\t\${num_diff_before}\t\${num_diff_after}"
      echo -e "${id}\tfilter_germline\text\t\${num_ext_before}\t\${num_ext_after}"
    } > ${id}.filter_germline.regions.tsv


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
      echo -e "${id}\tfilter_germline\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.filter_germline.metrics.tsv



    """
}

