process split_vcf {
    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi)

    output:
    tuple val(id), path("chunk_*.vcf.gz"), path("chunk_*.vcf.gz.tbi"), emit: vcf

    script:
    """
    # Extract contigs
    bcftools index -s ${vcf} > contigs.list

    ## Number of chunks to create
    #N=20

    ## Split contigs evenly into N groups
    #split -e -n l/\$N contigs.list contigs_

    # Create N sub-VCFs
    #for file in contigs_* ; do
    while read -r c p num; do
        chunk=\$c
        out=chunk_\$c.vcf.gz

        regions=\$c

        bcftools view -Oz -o \$out -r \$regions ${vcf}
        bcftools index -t \$out
    done < contigs.list
    """
}

