process run_phasing {

    publishDir "${params.results_dir}/phasing/intermediate_files",
    pattern: "${id}*.tsv*",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/",
    pattern: "${id}.phased.vcf.gz*",
    mode:'copy'

    cpus 8 
    memory '24G'
    time '6h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), 
        path(truth_vcf, stageAs: "?/*"), path(truth_vcf_tbi, stageAs: "?/*"),
        path(sr_bams), path(sr_bais), 
        path(lr_bams), path(lr_bais), 
        path(lr_ont_bams), path(lr_ont_bais),
        val(sex),
        path(germline_vcf), path(germline_tbi)
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)
    tuple path(ref), path(ref_index), path(ref_dict)


    output:
    tuple val(id), path("${id}.phased.vcf.gz"), path("${id}.phased.vcf.gz.tbi"),
        path(truth_vcf), path(truth_vcf_tbi), emit: vcf
    path("${id}.phasing.metrics.tsv")
    path("${id}.phasing.regions.tsv")
    path("${id}_phasing_tags.tsv.gz")
    path("${id}.read_counts.phasing.tsv")
    path("${id}.germline_map.tsv")

    script:
    """
    # Build: --pb-cram <bam1> --pb-cram <bam2> ...
    pb_crams=""
    for f in ${lr_bams}; do
        pb_crams+=" --pb-cram \${f}"
    done


    phase_mosaic_vars.sh -i ${id} \
        -v ${germline_vcf} \
        -w ${vcf} \
        \${pb_crams} \
        -r ${ref} -s ${sex} -t ${task.cpus}


    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.phased.vcf.gz
    # check regions
    num_easy_before=\$( bedtools intersect -u -b $easy_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_easy_after=\$( bedtools intersect -u -b $easy_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_diff_before=\$( bedtools intersect -u -b $diff_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_diff_after=\$( bedtools intersect -u -b $diff_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_ext_before=\$( bedtools intersect -u -b $ext_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_ext_after=\$( bedtools intersect -u -b $ext_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    {
      echo -e "id\tstep\tregion_type\tnum_before\tnum_after"
      echo -e "${id}\tphasing\teasy\t\${num_easy_before}\t\${num_easy_after}"
      echo -e "${id}\tphasing\tdiff\t\${num_diff_before}\t\${num_diff_after}"
      echo -e "${id}\tphasing\text\t\${num_ext_before}\t\${num_ext_after}"
    } > ${id}.phasing.regions.tsv


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
      echo -e "${id}\tphasing\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.phasing.metrics.tsv


    """
}
