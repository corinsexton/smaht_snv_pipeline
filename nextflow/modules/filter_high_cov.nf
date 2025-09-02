process filter_high_cov {

    publishDir "${params.results_dir}/high_cov_filtered",
    pattern: "${id}.filtered.vcf.gz*",
    mode:'copy'


    cpus 1
    memory '2G'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), path(bam), path(bai)
    val max_depth

    output:
    tuple val(id),
          path("${id}.filtered.vcf.gz"),
          path("${id}.filtered.vcf.gz.tbi")

    script:
    """
    filter_cov_lt_2x_mean.py --baseq 30 --threshold ${max_depth} ${vcf} ${bam} ${id}.filtered.vcf.gz
    tabix -f ${id}.filtered.vcf.gz
    """
}

