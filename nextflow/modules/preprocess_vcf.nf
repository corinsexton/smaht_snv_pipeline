
// norm, filter by PASS, atomize, and remove any duplicates from SNV caller vcf files.

process preprocess_vcf {
    cpus 1
    memory '4G'
    time '30m'

    publishDir "${params.results_dir}/1_pass_filtered",
    pattern: "${id}.preprocess.*.tsv",
    mode:'copy'

    input:
    tuple val(id), val(caller), path(vcf), path(tbi), path(ref), path(ref_index), path(truth_vcf), path(truth_tbi)
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)

    output:
    tuple val(id), val(caller), path("${id}.${caller}.norm.PASS.atom.dedup.vcf.gz"), path("${id}.${caller}.norm.PASS.atom.dedup.vcf.gz.tbi"), path(truth_vcf), path(truth_tbi), emit: vcf
    path("${id}.${caller}.preprocess.metrics.tsv"), emit: metrics
    path("${id}.preprocess.regions.tsv"), emit: regions


    script:
    """

    bcftools_PASS_norm_dedup.sh -i ${vcf} -r ${ref} -o ${id}.${caller}.norm.PASS.atom.dedup
    
    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.${caller}.norm.PASS.atom.dedup.vcf.gz

    num_before=\$(bcftools view -H "\${BEFORE_VCF}" | wc -l | awk '{print \$1}')
    num_after=\$(bcftools view -H "\${AFTER_VCF}"  | wc -l | awk '{print \$1}')

    ## check regions
    #num_easy_before=\$( bedtools intersect -u -b $easy_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    #num_easy_after=\$( bedtools intersect -u -b $easy_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    #num_diff_before=\$( bedtools intersect -u -b $diff_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    #num_diff_after=\$( bedtools intersect -u -b $diff_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    #num_ext_before=\$( bedtools intersect -u -b $ext_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    #num_ext_after=\$( bedtools intersect -u -b $ext_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    #{
    #  echo -e "id\tstep\tregion_type\tnum_before\tnum_after"
    #  echo -e "${id}\tpreprocess\teasy\t\${num_easy_before}\t\${num_easy_after}"
    #  echo -e "${id}\tpreprocess\tdiff\t\${num_diff_before}\t\${num_diff_after}"
    #  echo -e "${id}\tpreprocess\text\t\${num_ext_before}\t\${num_ext_after}"
    #} > ${id}.preprocess.regions.tsv
    touch ${id}.preprocess.regions.tsv


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
      echo -e "${id}\tpreprocess_${caller}\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.${caller}.preprocess.metrics.tsv


    """
}

