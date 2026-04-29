process run_minipileup2_parallel {

    cache 'lenient'

    cpus 1
    memory '24G'
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

    # Build PB long-read arrays for paired --lr-cram/--lr-tissue/--lr-type args
    pb_cram_arr=()
    for f in ${lr_bams}; do
        pb_cram_arr+=("\${f}")
    done

    lr_tissues_clean=\$(printf '%s\n' "${lr_tissues}" | tr -d '[],')
    pb_tissue_arr=()
    for f in \${lr_tissues_clean}; do
        pb_tissue_arr+=("\${f}")
    done

    pb_lr_args=""
    for i in "\${!pb_cram_arr[@]}"; do
        pb_lr_args+=" --lr-cram \${pb_cram_arr[\${i}]} --lr-tissue \${pb_tissue_arr[\${i}]} --lr-type PB"
    done

    # Build ONT long-read arrays for paired --lr-cram/--lr-tissue/--lr-type args
    ont_cram_arr=()
    for f in ${lr_ont_bams}; do
        ont_cram_arr+=("\${f}")
    done

    ont_tissues_clean=\$(printf '%s\n' "${lr_ont_tissues}" | tr -d '[],')
    ont_tissue_arr=()
    for f in \${ont_tissues_clean}; do
        ont_tissue_arr+=("\${f}")
    done

    ont_lr_args=""
    for i in "\${!ont_cram_arr[@]}"; do
        ont_lr_args+=" --lr-cram \${ont_cram_arr[\${i}]} --lr-tissue \${ont_tissue_arr[\${i}]} --lr-type ONT"
    done

    chr=\$( basename -s .vcf.gz ${vcf})

    minipileup2-parallel.sh -i ${vcf} \
        -r ${ref} \
        -o ${id}.\${chr}.minipileup \
        \${sr_crams} \
        \${pb_lr_args} \
        \${ont_lr_args}

    """
}
