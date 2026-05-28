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

    joined = vcf_inputs.join(bam_inputs)

    // Samples without PacBio LR data pass through unphased
    no_lr_passthrough = joined
       .filter { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_tissues, lr_ont_bams, lr_ont_bais, lr_ont_tissues -> !lr_bams }
       .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_tissues, lr_ont_bams, lr_ont_bais, lr_ont_tissues ->
           tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi)
       }

    joined
       .filter { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_tissues, lr_ont_bams, lr_ont_bais, lr_ont_tissues -> lr_bams as boolean }
       .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_tissues, lr_ont_bams, lr_ont_bais, lr_ont_tissues ->
         tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais)
       }.join(sex_ch)
       .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais, sex ->
         tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais, sex )
       }.join(germline_inputs)
       .set { phasing_input }

    run_phasing(phasing_input,regions_input,ref_input)

    all_output = run_phasing.out.vcf.mix(no_lr_passthrough)

    emit:
    all_output
}
