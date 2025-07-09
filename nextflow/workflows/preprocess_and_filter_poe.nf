/*
 * Workflow to preprocess VCFs and filter using panel of errors (POE)
 * 
 * Steps:
 *  1. Normalize, filter PASS, and atomize VCF (preprocess_vcf)
 *  2. Filter out known artifact variants from POE VCF (filter_panel_errors)
 */

nextflow.enable.dsl=2

include { preprocess_vcf }       from '../modules/preprocess_vcf.nf'
include { filter_panel_errors }  from '../modules/filter_panel_errors.nf'

workflow preprocess_and_filter_poe {


    take:
        vcf_inputs
        panel_of_errors_vcf
        ref

    main: 
    // Step 1: Normalize, filter PASS, atomize
    // Assumes: preprocess_vcf emits (id, processed_vcf, processed_vcf.tbi)
    vcf_inputs
        .map { id, vcf, tbi ->
            tuple(id, vcf, tbi, ref)
        }
        | preprocess_vcf

    // Step 2: Add POE + its index to each tuple
    filter_panel_errors(preprocess_vcf.out,panel_of_errors_vcf.first())

    emit:
    filter_panel_errors.out
}
