process phasing_step5 {

    publishDir "${params.results_dir}/phasing/step5/intermediate",
    pattern: "${id}.*.tsv",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/step5/failed",
    pattern: "${id}.phased.fail.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/step5/passed",
    pattern: "${id}.phased.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/phasing/step5/metrics",
    pattern: "${id}.phasing.metrics.tsv",
    mode:'copy'

    cpus 1
    memory '2G'
    time '1h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), 
        path(truth_vcf), path(truth_vcf_tbi),
        path(sr_bams), path(sr_bais), 
        path(lr_bams), path(lr_bais), 
        path(lr_ont_bams), path(lr_ont_bais),
        path(step4_tsv)


    output:
    tuple val(id), path("${id}.phased.vcf.gz"), path("${id}.phased.vcf.gz.tbi"),
        path(truth_vcf), path(truth_vcf_tbi), emit: vcf
    path("${id}.phasing.metrics.tsv"), emit: metrics

    script:
    """
    phasing_step5.py -t ${step4_tsv} -b ${lr_bams} -s 'male' -i ${id}

    sort -k1,1V -k2,2n ${id}_phasing_tags.tsv | bgzip - -o ${id}_phasing_tags.tsv.gz
    tabix -s1 -b2 -e2 ${id}_phasing_tags.tsv.gz

    bcftools annotate -a ${id}_phasing_tags.tsv.gz -c CHROM,POS,PHASING \
         -H '##INFO=<ID=PHASING,Number=1,Type=String,Description="Phasing classification from long-read haplotyping">' \
         -Oz -o annotated.vcf.gz ${vcf}
    tabix annotated.vcf.gz

    bcftools view -i '(FILTER="TIER1") || (INFO/PHASING="MOSAIC_PHASED") || (INFO/PHASING="UNABLE_TO_PHASE")' annotated.vcf.gz -Oz -o ${id}.phased.vcf.gz
    tabix ${id}.phased.vcf.gz

    bcftools view -i '(INFO/PHASING="ARTIFACT") || (INFO/PHASING="GERMLINE")' annotated.vcf.gz -Oz -o ${id}.phased.fail.vcf.gz
    tabix ${id}.phased.fail.vcf.gz


    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.phased.vcf.gz

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
