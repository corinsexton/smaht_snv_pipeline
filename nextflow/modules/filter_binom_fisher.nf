process filter_binom_fisher {

    publishDir "${params.results_dir}/binom_fisher",
    pattern: "${id}.tiered.binom_fisher.vcf.gz*",
    mode:'copy'

    cpus 1
    memory '1G'
    time '1h'

    tag "$id"

    input:
    tuple val(id), path(tiered_vcf), path(tiered_vcf_index)

    output:
    tuple val(id), path("${id}.tiered.binom_fisher.vcf.gz"),path("${id}.tiered.binom_fisher.vcf.gz.tbi")

    script:
    """
    filter_binom_fisher.py ${tiered_vcf} ${id}.tiered.binom_fisher.vcf
    bgzip ${id}.tiered.binom_fisher.vcf
    tabix ${id}.tiered.binom_fisher.vcf.gz
    """
}
