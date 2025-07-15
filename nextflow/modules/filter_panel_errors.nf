process filter_panel_errors {

    publishDir "${params.results_dir}/poe_filtered",
    pattern: "${id}.filtered.vcf.gz*",
    mode:'copy'


    cpus 1
    memory '2G'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi)
    tuple path(error_panel_vcf),path(error_panel_tbi)   // input index for panel
    path ucsc_regions
    path centromere_regions

    output:
    tuple val(id),
          path("${id}.filtered.vcf.gz"),
          path("${id}.filtered.vcf.gz.tbi")

    script:
    """
    filter_panel_errors.py ${vcf} ${error_panel_vcf} ${id}.filtered.vcf.gz

    # Exclude variants in UCSC SegDup regions
    bcftools view -T ^${ucsc_regions} ${id}.filtered.vcf.gz | \
        bcftools view -T ^${centromere_regions} -Oz -o filtered.vcf.gz

    mv filtered.vcf.gz ${id}.filtered.vcf.gz
    tabix -f ${id}.filtered.vcf.gz
    """
}

