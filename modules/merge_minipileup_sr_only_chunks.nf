process merge_minipileup_sr_only_chunks {

    publishDir "${params.results_dir}/12_minipileup_sr",
    pattern: "${id}.minipileup_sr.merged.vcf.gz*",
    mode: 'copy'

    cache 'lenient'

    cpus 1
    memory '2G'
    time '30m'

    tag "$id"

    input:
    tuple val(id), path(mp_chunk_vcfs), path(mp_chunk_tbis),
          path(orig_vcf), path(orig_tbi),
          path(truth_vcf), path(truth_tbi)

    output:
    tuple val(id),
          path(orig_vcf), path(orig_tbi),
          path(truth_vcf), path(truth_tbi),
          path("${id}.minipileup_sr.merged.vcf.gz"),
          path("${id}.minipileup_sr.merged.vcf.gz.tbi")

    script:
    """
    # Sort chunk files by name to ensure chromosome order
    readarray -t sorted_chunks < <(printf '%s\n' ${mp_chunk_vcfs.join(' ')} | sort -V)

    bcftools concat -Oz -o ${id}.minipileup_sr.merged.vcf.gz \\
        "\${sorted_chunks[@]}"
    bcftools index -t ${id}.minipileup_sr.merged.vcf.gz
    """
}
