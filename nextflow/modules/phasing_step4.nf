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
        path(truth_vcf, stageAs: "?/*"), path(truth_vcf_tbi, stageAs: "?/*"),
        path(germline_vcf), path(germline_tbi)

    output:
    tuple val(id), path("${id}.germline_map.tsv")

    script:
    """
    phasing_step1_get_closest_germline.sh -s $id -v $germline_vcf -w $vcf
    """
}
