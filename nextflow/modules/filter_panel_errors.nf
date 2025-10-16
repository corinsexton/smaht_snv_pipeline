process filter_poe {

    publishDir "${params.results_dir}/poe_filtered",
    pattern: "${id}.pon.filtered.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/poe_filtered",
    pattern: "${id}.filter_poe.*.tsv",
    mode:'copy'

    cache 'lenient'
    cpus 1
    memory '8G'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi)
    tuple path(error_panel_fa),path(error_panel_fai)   // input index for panel
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)

    output:
    tuple val(id),
          path("${id}.pon.filtered.vcf.gz"),
          path("${id}.pon.filtered.vcf.gz.tbi"),
          path(truth_vcf), path(truth_tbi), emit: vcf
    path("${id}.filter_poe.metrics.tsv"), emit: metrics
    path("${id}.filter_poe.regions.tsv"), emit: regions

    script:
    """
    filter_by_poe.py --vcf ${vcf} \
                       --fasta ${error_panel_fa} \
                       --out ${id}.pon.filtered.vcf 
    # optional params
    #--threads 2 --failed-out failed.vcf

    bcftools view -Oz ${id}.pon.filtered.vcf > ${id}.pon.filtered.vcf.gz
    tabix ${id}.pon.filtered.vcf.gz

     # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.pon.filtered.vcf.gz

    # check regions
    num_easy_before=\$( bedtools intersect -u -b $easy_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_easy_after=\$( bedtools intersect -u -b $easy_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_diff_before=\$( bedtools intersect -u -b $diff_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_diff_after=\$( bedtools intersect -u -b $diff_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_ext_before=\$( bedtools intersect -u -b $ext_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_ext_after=\$( bedtools intersect -u -b $ext_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')


    {
      echo -e "id\tstep\tregion_type\tnum_before\tnum_after"
      echo -e "${id}\tfilter_poe\teasy\t\${num_easy_before}\t\${num_easy_after}"
      echo -e "${id}\tfilter_poe\tdiff\t\${num_diff_before}\t\${num_diff_after}"
      echo -e "${id}\tfilter_poe\text\t\${num_ext_before}\t\${num_ext_after}"
    } > ${id}.filter_poe.regions.tsv



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
      echo -e "${id}\tfilter_poe\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.filter_poe.metrics.tsv




    """
}

