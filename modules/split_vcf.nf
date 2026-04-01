process split_vcf {
    tag "$id"

    cpus 1
    memory '1G'
    time '30m'

    input:
    tuple val(id), path(vcf), path(tbi)
    path(chunks_file)

    output:
    tuple val(id), path("chunk_*.vcf.gz"), path("chunk_*.vcf.gz.tbi"), emit: vcf

    script:
    """
    idx=0
    while IFS= read -r region || [[ -n "\$region" ]]; do
        [[ -z "\$region" ]] && continue
        idx=\$((idx+1))
        out=\$(printf "chunk_%03d.vcf.gz" \$idx)
        bcftools view -Oz -o "\$out" -r "\$region" ${vcf}
        bcftools index -t "\$out"
    done < "${chunks_file}"
    """
}
