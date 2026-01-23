process merge_minipileup_chunks {
    tag "$id"

    input:
    tuple val(id), path(mp_chunk_vcfs), path(mp_chunk_tbis),
          path(orig_vcf), path(orig_tbi),
          path(truth_vcf), path(truth_tbi)

    output:
    tuple val(id),
          path(orig_vcf), path(orig_tbi),
          path(truth_vcf), path(truth_tbi),
          path("${id}.minipileup.merged.vcf.gz")

    script:
    """
    bcftools concat -Oz -o ${id}.minipileup.merged.vcf.gz \\
        ${chunk_vcfs.join(' ')}
    bcftools index -t ${id}.minipileup.merged.vcf.gz
    """
}

