nextflow.enable.dsl=2

include { run_phasing } from '../modules/run_phasing.nf'

workflow phasing {

    take:
        vcf_inputs
        germline_inputs
        bam_inputs
        ref_input
        vep_config
        regions_input
        sex_ch

    main: 

    vcf_inputs
       .join(bam_inputs)     // join on 'id'
       .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_tissues, lr_ont_bams, lr_ont_bais, lr_ont_tissues ->
         tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais)
       }.join(sex_ch)
       .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais, sex ->
         tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais, sex )
       }.join(germline_inputs)
       .set { phasing_input }

    run_phasing(phasing_input,regions_input,ref_input)

    emit:
    run_phasing.out.vcf
}
