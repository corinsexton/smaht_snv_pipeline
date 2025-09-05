process tier_variants {

    publishDir "${params.results_dir}/tiered_variants",
    pattern: "${id}.tiered.vcf.gz*",
    mode:'copy'

    publishDir "${params.results_dir}/minipileup/parsed",
    pattern: "${id}.parsed.minipileup.tsv",
    mode:'copy'

    cpus 1
    memory '1G'
    time '1h'

    tag "$id"

    input:
    tuple val(id), 
          path(vcf), path(tbi),
          path(minipileup_vcf)

    output:
    tuple val(id), path("${id}.tiered.vcf.gz"), path("${id}.tiered.vcf.gz.tbi")

    script:
    """
    # get the counts for pileups
    # matches ref and alt alleles
    parse_minipileup.py ${vcf} \
            ${minipileup_vcf} \
            ${id}.parsed.minipileup.tsv

    # assigns tiers based on LR ≥ 1 (TIER1), LR = 0 & SR ≥ 2 (TIER2), else (.)
    split_lr_presence.py ${id}.parsed.minipileup.tsv ${id}.tiered.vcf
    bgzip ${id}.tiered.vcf
    tabix ${id}.tiered.vcf.gz
    """
}
