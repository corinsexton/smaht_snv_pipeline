process phasing_step4 {

    publishDir "${params.results_dir}/phasing/step4",
    pattern: "${id}.germline_map.tsv",
    mode:'copy'

    cpus 1
    memory '2G'
    time '1h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi),
        path(truth_vcf), path(truth_vcf_tbi)
    path(windows_bed)

    output:
    tuple val(id), path("${id}.germline_map.tsv")

    script:
    """
    phasing_step4.sh -s $id -v $vcf -w $windows_bed
    """
}
