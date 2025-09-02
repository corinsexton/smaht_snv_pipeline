process filter_centromere_segdups {

    publishDir "${params.results_dir}/centromere_segdups_filtered",
    pattern: "${id}.filtered.vcf.gz*",
    mode:'copy'


    cpus 1
    memory '2G'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi)
    path ucsc_regions
    path centromere_regions

    output:
    tuple val(id),
          path("${id}.filtered.vcf.gz"),
          path("${id}.filtered.vcf.gz.tbi")

    script:
    """
    # Exclude variants in UCSC SegDup regions, centromeres, and >2 x avg_depth
    bcftools view -T ^${ucsc_regions} ${vcf} | \
        bcftools view -T ^${centromere_regions} -Oz -o filtered.vcf.gz

    mv filtered.vcf.gz ${id}.filtered.vcf.gz
    tabix -f ${id}.filtered.vcf.gz
    """
}

