nextflow.enable.dsl=2

include { run_minipileup2_sr_only_parallel } from '../modules/run_minipileup2_sr_only_parallel.nf'
include { run_minipileup_sr_only }           from '../modules/run_minipileup_sr_only.nf'

workflow check_other_tissues {

    take:
        vcf_inputs
        ref_input
        bam_inputs
        regions_input

    main:

    vcf_inputs
      .join(bam_inputs)
      .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, tissue_ids ->
        tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, tissue_ids)
      }
      .set { vcf_bam_channel }

    run_minipileup2_sr_only_parallel(vcf_bam_channel, ref_input, regions_input)

    run_minipileup_sr_only(run_minipileup2_sr_only_parallel.out.vcf, ref_input, regions_input)

    emit:
    run_minipileup_sr_only.out.vcf
}
