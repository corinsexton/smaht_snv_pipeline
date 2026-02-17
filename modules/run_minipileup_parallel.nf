process run_minipileup_parallel {

    cache 'lenient'

    cpus 10
    memory '4G'
    time '2h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), 
        path(truth_vcf), path(truth_vcf_tbi),
        path(sr_bams), path(sr_bais), 
        path(lr_bams), path(lr_bais), val(lr_tissues),
        path(lr_ont_bams), path(lr_ont_bais), val(lr_ont_tissues)
    tuple path(ref), path(ref_index), path(ref_dict)


    output:
    tuple val(id),
          path(vcf), path(tbi), path(truth_vcf), path(truth_vcf_tbi),
          path("${id}*.minipileup.vcf.gz"), path("${id}*.minipileup.vcf.gz.tbi"), emit: vcf

    script:
    """

    # Build: --sr-cram <bam1> --sr-cram <bam2> ...
    sr_crams=""
    for f in ${sr_bams}; do
        sr_crams+=" --sr-cram \${f}"
    done

    # Build: --pb-cram <bam1> --pb-cram <bam2> ...
    pb_crams=""
    for f in ${lr_bams}; do
        pb_crams+=" --pb-cram \${f}"
    done

    # Build: --ont-cram <bam1> --ont-cram <bam2> ...
    ont_crams=""
    for f in ${lr_ont_bams}; do
        ont_crams+=" --ont-cram \${f}"
    done

    lr_tissues_clean=\$(printf '%s\n' "${lr_tissues}" | tr -d '[],')
    # Build: --pb-tissue <tissue1> --pb-tissue <tissue2> ...
    pb_tissues=""
    for f in \${lr_tissues_clean}; do
        pb_tissues+=" --pb-tissue \${f}"
    done

    ont_tissues_clean=\$(printf '%s\n' "${lr_ont_tissues}" | tr -d '[],')
    # Build: --ont-tissue <tissue1> --ont-tissue <tissue2> ...
    ont_tissues=""
    for f in \${ont_tissues_clean}; do
        ont_tissues+=" --ont-tissue \${f}"
    done

    chr=\$( basename -s .vcf.gz ${vcf})

    minipileup-parallel.sh -i ${vcf} \
        -r ${ref} \
        -t ${task.cpus} \
        --group 50 \
        -o ${id}.\${chr}.minipileup \
        \${sr_crams} \
        \${pb_crams} \
        \${ont_crams} \
        \${pb_tissues} \
        \${ont_tissues}

    """
}
