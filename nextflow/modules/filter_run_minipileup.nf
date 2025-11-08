process filter_run_minipileup {

    publishDir "${params.results_dir}/minipileup",
    pattern: "${id}.minipileup.vcf",
    mode:'copy'

    cache 'lenient'


    cpus 20
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

    export ont_bams
    # --- Function to prepare region-specific BAM inputs (only convert CRAMs) ---
    build_bam_inputs() {
        local chr="\$1"
        local start="\$2"
        local end="\$3"
        shift 3
        local bam_list=()
        mkdir -p /n/data1/hms/dbmi/park-smaht_dac/analysis/tmp/${id}
        for f in "\$@"; do
            if [[ "\$f" == *.cram ]]; then
                tmp_bam="\$(mktemp -p /n/data1/hms/dbmi/park-smaht_dac/analysis/tmp/${id} --suffix=.bam)"
                samtools view -T ${ref} -b "\$f" "\${chr}:\${end}-\${end}" -o "\$tmp_bam"
                samtools index "\$tmp_bam"
                bam_list+=("\$tmp_bam")
            else
                bam_list+=("\$f")
            fi
        done
        echo "\${bam_list[*]}"
    }

    export -f build_bam_inputs



    # First extract intervals
    bcftools query -f '%CHROM\t%POS0\t%END\n' "${vcf}" > "${id}.bed"


    read -r chr start end < "${id}.bed"
    bam_inputs=\$(build_bam_inputs "\$chr" "\$end" "\$end" ${sr_bams} ${lr_bams} \${ont_bams})
    echo here
    echo \$bam_inputs
    minipileup -f "${ref}" -c -C -T 5 -Q 30 -q 10 -r "\${chr}:\${end}-\${end}" \${bam_inputs} > x    # weirdly not working when crams, hacky fix
    grep '^#' x > "${id}.minipileup.vcf"
    echo here2
    for bam in \${bam_inputs}; do
        if [[ \$bam == */tmp/*bam ]]; then
            rm "\$bam"
            rm "\${bam}.bai"
       fi
    done

    for bam in ${sr_bams} ${lr_bams} \${ont_bams}; do
        echo \$bam
    done > new_samples.txt

    bcftools query -l "${id}.minipileup.vcf" > old_samples.txt
    paste old_samples.txt new_samples.txt > rename.txt
    bcftools reheader -s rename.txt -o rename.vcf "${id}.minipileup.vcf"
    mv rename.vcf "${id}.minipileup.vcf"


    # Run each interval in parallel with 8 jobs
    # BQ ≥ 30, MQ ≥ 20, count alleles both strands (-C), vcf format (-c), trim 5bp each end (-T 5)
    # -s drop alleles with depth<INT (0)
    echo here3
    cat "${id}.bed" | xargs -P20 -n3 bash -c '
        chr=\$1; pos1=\$2; pos2=\$3
        tmp=\$(mktemp)
        tmp2=\$(mktemp)
        bam_inputs=\$(build_bam_inputs "\$chr" "\$pos2" "\$pos2" ${sr_bams} ${lr_bams} \${ont_bams})
        minipileup -f "'"${ref}"'" -c -C -Q 20 -q 30 -T 5 -s 0 -r "\${chr}:\${pos2}-\${pos2}" \${bam_inputs} > \${tmp2}
        grep -v "^#" \${tmp2} > "\$tmp"
        rm \${tmp2}
        for bam in \${bam_inputs}; do
            if [[ \$bam == */tmp/*bam ]]; then
                rm "\$bam"
                rm "\${bam}.bai"
            fi
        done
        echo "\$tmp"
    ' _ > tmp_files.list
    
    # concatenate, sort, dedup
    cat \$(cat tmp_files.list) >> "${id}.minipileup.vcf"
    
    # cleanup
    xargs rm -f < tmp_files.list
    rm -f tmp_files.list
    rm -rf /n/data1/hms/dbmi/park-smaht_dac/analysis/tmp/${id}

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
