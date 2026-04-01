process split_vcf {
    tag "$id"

    cpus 1
    memory '1G'
    time '30m'

    input:
    tuple val(id), path(vcf), path(tbi)

    output:
    tuple val(id), path("chunk_*.vcf.gz"), path("chunk_*.vcf.gz.tbi"), emit: vcf

    script:
    """
    CHROMS=(chr{1..22} chrX chrY)


    # Create N sub-VCFs
    #for file in contigs_* ; do
    for c in "\${CHROMS[@]}"; do
        chunk=\$c
        out=chunk_\$c.vcf.gz

        regions=\$c

        bcftools view -Oz -o \$out -r \$regions ${vcf}
        bcftools index -t \$out
    done
    """
}

