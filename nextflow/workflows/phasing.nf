nextflow.enable.dsl=2

include { phasing_step1 } from '../modules/phasing_step1.nf'
include { phasing_step4 } from '../modules/phasing_step4.nf'
include { phasing_step5 } from '../modules/phasing_step5.nf'
include { runVEP as runVEPonVCF } from '../modules/run_vep.nf'
include { filter_vep } from '../modules/filter_vep.nf'
include { runVEP } from '../modules/run_vep.nf'

workflow phasing {

    take:
        vcf_inputs
        germline_inputs
        bam_inputs
        ref_input
        vep_config
        regions_input

    main: 
    // Step 1: run minipileup to get counts in SR and LR
    vcf_inputs
      .join(bam_inputs)     // join on 'id'
      .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais ->
        tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais)
      }
      .set { vcf_bam_channel }

   //phasing_step1(vcf_bam_channel,ref_input)


   //// use VEP for AF filtering
   //phasing_step1.out.vcf.combine( vep_config )
   //                    .map {
   //                          id, vcf, tbi, truth_vcf, truth_tbi, config ->
   //                           [id: id, file: vcf, index: tbi, truth_vcf:truth_vcf, truth_tbi:truth_tbi, vep_config: config]
   //                         }
   //                    .set {vep_input}


   //// run and filter with vep
   //runVEP(vep_input,'germline')
   //filter_vep(runVEP.out,'germline')

   phasing_step4(vcf_inputs.join(germline_inputs))

   vcf_inputs
      .join(bam_inputs)     // join on 'id'
      .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais ->
        tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais)
      }.join(phasing_step4.out)
      .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais, step4_tsv ->
        tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais, step4_tsv)
      }
      .set { step5_input }


   phasing_step5(step5_input,regions_input)

    emit:
    phasing_step5.out.vcf
}
