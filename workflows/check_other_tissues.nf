nextflow.enable.dsl=2

include { run_minipileup_sr_only }       from '../modules/run_minipileup_sr_only.nf'

workflow check_other_tissues {

    take:
        vcf_inputs
        ref_input
        bam_inputs
        regions_input

    main: 
    // Step 1: run minipileup to get counts in SR and LR
    vcf_inputs
      .join(bam_inputs)     // join on 'id'
      .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, tissue_ids ->
        tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, tissue_ids)
      }
      .set { vcf_bam_channel }

    run_minipileup_sr_only(vcf_bam_channel, ref_input, regions_input)


    emit:
    run_minipileup_sr_only.out.vcf
}
