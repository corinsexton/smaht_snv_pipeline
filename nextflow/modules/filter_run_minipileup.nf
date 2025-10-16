process filter_run_minipileup {

    publishDir "${params.results_dir}/minipileup",
    pattern: "${id}.minipileup.vcf",
    mode:'copy'

    cache 'lenient'


    cpus 8
    memory '2G'
    time '6h'

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
          path("${id}.minipileup.vcf"), emit: vcf
    path("labels.txt"), emit: labels

    script:
    """

     # --- Optional ONT input handling ---
    # If ONT file is named "none.bam", set empty variable
    ont_bams=""
    if [[ ! "${lr_ont_bams}" == *"none.bam"* ]]; then
        ont_bams="${lr_ont_bams}"
    fi

    # convert to bedfile for minipileup


    # First extract intervals
    bcftools query -f '%CHROM\t%POS0\t%END\n' "${vcf}" > "${id}.bed"

    echo 'here'
    

    read -r chr start end < "${id}.bed"
    minipileup -f "${ref}" \
        -c -C -T 5 -Q 30 -q 10 \
        -r "\${chr}:\${start}-\${end}" \
        ${sr_bams} ${lr_bams} \${ont_bams} \
        | grep '^#' > "${id}.minipileup.vcf"
    
    # Run each interval in parallel with 8 jobs
    # BQ ≥ 30, MQ ≥ 10, count alleles both strands (-C), vcf format (-c), trim 5bp each end (-T 5)
    # -s drop alleles with depth<INT (0)
    cat "${id}.bed" | xargs -P8 -n3 bash -c '
        chr=\$1; pos1=\$2; pos2=\$3
        tmp=\$(mktemp)
        minipileup -f "'"${ref}"'" \
            -c -C -Q 20 -q 30 \
            -s 0 \
            -r "\${chr}:\${pos2}-\${pos2}" \
            '"${sr_bams} ${lr_bams} \${ont_bams}"' \
        | grep -v "^#" > "\$tmp"
        echo "\$tmp"
    ' _ > tmp_files.list
    
    # concatenate, sort, dedup
    cat \$(cat tmp_files.list) >> "${id}.minipileup.vcf"
    
    # cleanup
    xargs rm -f < tmp_files.list
    rm -f tmp_files.list

    bcftools sort ${id}.minipileup.vcf


    # build label string
    labels=""
    for _ in ${sr_bams}; do labels+="SR,"; done
    for _ in ${lr_bams}; do labels+="LR,"; done
    for _ in \${ont_bams}; do labels+="ONT,"; done

    # remove trailing comma and save to file
    echo \${labels%,} > labels.txt
    """
}
