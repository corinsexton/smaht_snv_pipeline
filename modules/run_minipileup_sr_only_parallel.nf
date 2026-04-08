process run_minipileup_sr_only_parallel {

    cache 'lenient'

    cpus 15
    memory '4G'
    time '2h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi),
        path(truth_vcf), path(truth_vcf_tbi),
        path(sr_bams), path(sr_bais), val(sr_ids)
    tuple path(ref), path(ref_index), path(ref_dict)
    tuple path(easy_regions), path(diff_regions), path(ext_regions),
        path(easy_regions_tbi), path(diff_regions_tbi), path(ext_regions_tbi)

    output:
    tuple val(id),
          path(vcf), path(tbi), path(truth_vcf), path(truth_vcf_tbi),
          path("${id}*.minipileup_sr.vcf.gz"), path("${id}*.minipileup_sr.vcf.gz.tbi"), emit: vcf

    script:
    """
    # Build: --sr-cram <bam1> --sr-cram <bam2> ...
    sr_crams=""
    for f in ${sr_bams}; do
        sr_crams+=" --sr-cram \${f}"
    done

    # Build: --sr-tissue <tissue id1> --sr-tissue <tissue id2> ...
    sr_tissue=""
    for f in ${sr_ids}; do
        sr_tissue+=" --sr-tissue \${f}"
    done

    chunk=\$( basename -s .vcf.gz ${vcf})

    # Filter chunk to easy regions before running minipileup
    bcftools view -R ${easy_regions} -Oz -o \${chunk}.easyonly.vcf.gz ${vcf}
    tabix \${chunk}.easyonly.vcf.gz

    minipileup-parallel_sr_only.sh -i \${chunk}.easyonly.vcf.gz \
        -r ${ref} \
        -t ${task.cpus} \
        --group 50 \
        -o ${id}.\${chunk}.minipileup_sr \
        \${sr_crams} \
        \${sr_tissue}
    """
}
