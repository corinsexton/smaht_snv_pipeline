process filter_clustered_variants {

    publishDir "${params.results_dir}/clusters_filtered",
    pattern: "${id}.filtered.clusters.vcf.gz*",
    mode:'copy'


    publishDir "${params.results_dir}/clusters_filtered",
    pattern: "${id}.filter_clustered.metrics.tsv",
    mode:'copy'


    cpus 1
    memory '500M'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi)

    output:
    tuple val(id),
          path("${id}.filtered.clusters.vcf.gz"),
          path("${id}.filtered.clusters.vcf.gz.tbi"),
          path(truth_vcf), path(truth_tbi), emit: vcf
        path("${id}.filter_clustered.metrics.tsv"), emit: metrics

    script:
    """


    filter_clustered_variants.py --window 100 ${vcf} ${id}.filtered.clusters.vcf
    bgzip ${id}.filtered.clusters.vcf
    tabix ${id}.filtered.clusters.vcf.gz


    # --- metrics (standard schema) ---
    BEFORE_VCF=${vcf}
    AFTER_VCF=${id}.filtered.clusters.vcf.gz

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
      echo -e "${id}\tfilter_clustered\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.filter_clustered.metrics.tsv



    """
}

