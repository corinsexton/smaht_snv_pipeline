process run_minipileup_parallel {

    publishDir "${params.results_dir}/minipileup",
    pattern: "${id}.minipileup.vcf.gz*",
    mode:'copy'

    cache 'lenient'


    cpus 30
    memory '8G'
    time '12h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), 
        path(truth_vcf), path(truth_vcf_tbi),
        path(sr_bams), path(sr_bais), 
        path(lr_bams), path(lr_bais), 
        path(lr_ont_bams), path(lr_ont_bais)
    tuple path(ref), path(ref_index), path(ref_dict)


    output:
    tuple val(id),
          path(vcf), path(tbi), path(truth_vcf), path(truth_vcf_tbi),
          path("${id}.minipileup.vcf.gz"), emit: vcf

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

    minipileup-parallel.sh -i ${vcf} \
        -r ${ref} \
        -t 20 \
        -o ${id}.minipileup \
        \${sr_crams} \
        \${pb_crams} \
        \${ont_crams} 

    """
}
