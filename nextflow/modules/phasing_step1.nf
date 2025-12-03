process phasing_step1 {

    publishDir "${params.results_dir}/phasing/step1",
    pattern: "${id}.hc.norm.vcf.gz*",
    mode:'copy'

    cpus 4
    memory '6G'
    time '2h'

    tag "$id"

    input:
    tuple val(id), path(vcf), path(tbi),
        path(truth_vcf, stageAs: "?/*"), path(truth_vcf_tbi, stageAs: "?/*"),
        path(sr_bams), path(sr_bais),
        path(lr_bams), path(lr_bais),
        path(lr_ont_bams), path(lr_ont_bais)
    tuple path(ref), path(ref_index), path(ref_dict)

    output:
    tuple val(id),
        path("${id}.hc.norm.vcf.gz"),
        path("${id}.hc.norm.vcf.gz.tbi"),
        path(truth_vcf), path(truth_vcf_tbi), emit: vcf
    path("${id}.tier1_windows.bed"), emit: window_bed
        

    script:
    """
    phasing_step1.sh -r $ref -s $id -v $vcf -b $sr_bams -t 4
    """
}
