process filter_run_minipileup {

    publishDir "${params.results_dir}/minipileup",
    pattern: "${id}.minipileup.vcf",
    mode:'copy'


    cpus 1
    memory '2G'
    time '6h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi), 
        path(bam), path(bai), 
        path(lr_bam), path(lr_bai)
    tuple path(ref), path(ref_index)


    output:
    tuple val(id),
          path(vcf), path(tbi),
          path("${id}.minipileup.vcf")

    script:
    """
    # convert to bedfile for minipileup

    bcftools query -f '%CHROM\t%POS0\t%END\n' ${vcf} > ${id}.bed

    read -r chr start end < ${id}.bed

    minipileup -f ${ref} \
        -c -C -T 5 -Q 30 -q 10 \
        -r \$chr:\$start-\$end \
        ${bam} ${lr_bam} | grep '^#'  > ${id}.minipileup.vcf

    # # BQ ≥ 30, MQ ≥ 10, count alleles both strands (-C), vcf format (-c), trim 5bp each end (-T 5)
    # while read -r chr pos1 pos2; do
    #     line=\$( minipileup -f ${ref} \
    #         -c -C -T 5 -Q 30 -q 10 \
    #         -r \$chr:\$pos2-\$pos2 \
    #         ${bam} ${lr_bam} )
    #     if \$( echo \$line | grep -q -v '^#'); then
    #        \$( echo \$line | grep -v '^#')  >> ${id}.minipileup.vcf
    #     else
    #         echo \$chr \$pos2 \$pos2 >> failed.txt
    #     fi

    # BQ ≥ 30, MQ ≥ 10, count alleles both strands (-C), vcf format (-c), trim 5bp each end (-T 5)
    # -s drop alleles with depth<INT (0)
    while read -r chr pos1 pos2; do
        minipileup -f ${ref} \
            -c -C -Q 20 -q 30 \
            -s 0 \
            -r \$chr:\$pos2-\$pos2 \
            ${bam} ${lr_bam} | { grep -v '^#' || true; } >> ${id}.minipileup.vcf
    done < ${id}.bed
    """
}
