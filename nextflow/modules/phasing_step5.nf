process phasing_step5 {

    publishDir "${params.results_dir}/phasing/step5/intermediate",
    pattern: "${id}_phasing_tags.tsv.gz",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/step5/intermediate",
    pattern: "${id}.read_counts.phasing.tsv",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/step5/failed",
    pattern: "${id}.phased.fail.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/step5/passed",
    pattern: "${id}.phased.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/step5/metrics",
    pattern: "${id}.phasing.*.tsv",
    mode:'copy'

    cpus 20
    memory '16G'
    time '12h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), 
        path(truth_vcf), path(truth_vcf_tbi),
        path(sr_bams), path(sr_bais), 
        path(lr_bams), path(lr_bais), 
        path(lr_ont_bams), path(lr_ont_bais),
        path(step4_tsv), val(sex)
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)


    output:
    tuple val(id), path("${id}.phased.vcf.gz"), path("${id}.phased.vcf.gz.tbi"),
        path(truth_vcf), path(truth_vcf_tbi), emit: vcf
    path("${id}.phasing.metrics.tsv")
    path("${id}.phasing.regions.tsv")
    path("${id}_phasing_tags.tsv.gz")
    path("${id}.phased.fail.vcf.gz")
    path("${id}.read_counts.phasing.tsv")

    script:
    """
    phasing_step5.py -w 8 -t ${step4_tsv} -b ${lr_bams} -s ${sex} -i ${id}

    sort -k1,1V -k2,2n ${id}_phasing_tags.tsv | bgzip - -o ${id}_phasing_tags.tsv.gz
    tabix -s1 -b2 -e2 ${id}_phasing_tags.tsv.gz

    bcftools annotate -a ${id}_phasing_tags.tsv.gz -c CHROM,POS,PHASING \
         -H '##INFO=<ID=PHASING,Number=1,Type=String,Description="Phasing classification from long-read haplotyping">' \
         -Oz -o annotated.vcf.gz ${vcf}
    tabix annotated.vcf.gz

    bcftools view -i '(FILTER="TIER2") || (INFO/PHASING="MOSAIC_PHASED") || (INFO/PHASING="UNABLE_TO_PHASE")' annotated.vcf.gz -Oz -o ${id}.phased.vcf.gz
    tabix ${id}.phased.vcf.gz

    bcftools view -i '(INFO/PHASING="ARTIFACT") || (INFO/PHASING="GERMLINE")' annotated.vcf.gz -Oz -o ${id}.phased.fail.vcf.gz
    tabix ${id}.phased.fail.vcf.gz


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
