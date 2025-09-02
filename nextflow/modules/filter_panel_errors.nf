process filter_poe {

    publishDir "${params.results_dir}/poe_filtered",
    pattern: "${id}.pon.filtered.vcf.gz*",
    mode:'copy'


    cpus 1
    memory '500M'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi)
    tuple path(error_panel_fa),path(error_panel_fai)   // input index for panel

    output:
    tuple val(id),
          path("${id}.pon.filtered.vcf.gz"),
          path("${id}.pon.filtered.vcf.gz.tbi")

    script:
    """
    filter_by_poe.py --vcf ${vcf} \
                       --fasta ${error_panel_fa} \
                       --out ${id}.pon.filtered.vcf 
    # optional params
    #--threads 2 --failed-out failed.vcf

    bcftools view -Ob ${id}.pon.filtered.vcf > ${id}.pon.filtered.vcf.gz
    bcftools index ${id}.pon.filtered.vcf.gz
    """
}

