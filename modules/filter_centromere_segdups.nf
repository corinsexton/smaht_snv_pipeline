process filter_centromere_segdups {

    publishDir "${params.results_dir}/5_centromere_segdups_filtered",
    pattern: "${id}.filtered.vcf.gz*",
    mode:'copy'


    publishDir "${params.results_dir}/5_centromere_segdups_filtered",
    pattern: "${id}.filter_centromere_segdups.*.tsv",
    mode:'copy'


    cpus 1
    memory '4G'
    time '60m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi)
    path ucsc_regions
    path centromere_regions
    path simple_repeat_regions 
    path kg_indels
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)

    output:
        tuple val(id),
          path("${id}.filtered.vcf.gz"),
          path("${id}.filtered.vcf.gz.tbi"),
          path(truth_vcf), path(truth_tbi), emit: vcf
        path("${id}.filter_centromere_segdups.metrics.tsv"), emit: metrics
        path("${id}.filter_centromere_segdups.regions.tsv"), emit: regions 

    script:
    """

    filter_regions_keep_fails.py --segdup ${ucsc_regions} \
                                 --centromere ${centromere_regions} \
                                 --simple-repeat ${simple_repeat_regions} \
                                 --kg-indels ${kg_indels} --kg-slop 5 \
                                 ${vcf} ${id}.filtered.vcf.gz 

    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.filtered.vcf.gz

    # check regions
    num_easy_before=\$( bedtools intersect -u -b $easy_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_easy_after=\$( bedtools intersect -u -b $easy_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_diff_before=\$( bedtools intersect -u -b $diff_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_diff_after=\$( bedtools intersect -u -b $diff_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')

    num_ext_before=\$( bedtools intersect -u -b $ext_regions -a \${BEFORE_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')
    num_ext_after=\$( bedtools intersect -u -b $ext_regions -a \${AFTER_VCF} | grep -v "^#" | wc -l | awk '{print \$1}')


    {
      echo -e "id\tstep\tregion_type\tnum_before\tnum_after"
      echo -e "${id}\tfilter_centromere_segdups\teasy\t\${num_easy_before}\t\${num_easy_after}"
      echo -e "${id}\tfilter_centromere_segdups\tdiff\t\${num_diff_before}\t\${num_diff_after}"
      echo -e "${id}\tfilter_centromere_segdups\text\t\${num_ext_before}\t\${num_ext_after}"
    } > ${id}.filter_centromere_segdups.regions.tsv



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
      echo -e "${id}\tfilter_centromere_segdups\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.filter_centromere_segdups.metrics.tsv



    """
}

