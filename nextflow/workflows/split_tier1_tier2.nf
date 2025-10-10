nextflow.enable.dsl=2

include { filter_run_minipileup }       from '../modules/filter_run_minipileup.nf'
include { tier_variants }       from '../modules/tier_variants.nf'
include { filter_binom_fisher }       from '../modules/filter_binom_fisher.nf'

workflow split_tier1_tier2 {

    take:
        vcf_inputs
        bam_inputs
        ref_input
        regions_input

    main: 
    // Step 1: run minipileup to get counts in SR and LR
    vcf_inputs
      .join(bam_inputs)     // join on 'id'
      .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais ->
        tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais)
      }
      .set { vcf_bam_channel }

    filter_run_minipileup(vcf_bam_channel,ref_input)

    // count up reads in each assay, tier1/2 label
    tier_variants(filter_run_minipileup.out.vcf,filter_run_minipileup.out.labels,regions_input)

    // run binomial and fisher tests to get passing variants
    filter_binom_fisher(tier_variants.out.vcf, regions_input)

    emit:
    filter_binom_fisher.out.vcf
}
