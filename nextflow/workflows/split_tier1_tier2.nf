nextflow.enable.dsl=2

include { filter_run_minipileup }       from '../modules/filter_run_minipileup.nf'

workflow split_tier1_tier2 {

    take:
        vcf_inputs
        bam_inputs
        ref_input

    main: 
    // Step 1: run minipileup to get counts in SR and LR
    vcf_inputs
      .join(bam_inputs)     // join on 'id'
      .map { id, vcf, tbi, bam, bai, lr_bam, lr_bai ->
        tuple(id, vcf, tbi, bam, bai, lr_bam, lr_bai)
      }
      .set { vcf_bam_channel }

    filter_run_minipileup(vcf_bam_channel,ref_input)

    //    // input: tuple val(id), path(vcf), path(tbi), path("${id}.minipileup.vcf")
    //    // output: tuple val(id), path(vcf), path(tbi), path("${id}.minipileup.vcf")
    //    filter_depth(filter_run_minipileup.out) 
    //    filter_strand_bias(filter_depth.out) 


    //    // output: filter_tier1_tier2.tier1 tuple val(id),  path(vcf), path(tbi)
    //    // output: filter_tier1_tier2.tier2 tuple val(id),  path(vcf), path(tbi)
    //    filter_tier1_tier2(filter_strand_bias.out)

    emit:
    filter_run_minipileup.out
}
