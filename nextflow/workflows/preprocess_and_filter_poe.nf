/*
 * Workflow to preprocess VCFs and filter using panel of errors (POE)
 *  Filter out known artifact variants from POE VCF (filter_panel_errors)
 */

nextflow.enable.dsl=2

include { preprocess_vcf }       from '../modules/preprocess_vcf.nf'
include { filter_centromere_segdups }  from '../modules/filter_centromere_segdups.nf'
include { filter_poe }  from '../modules/filter_panel_errors.nf'
include { filter_clustered_variants }  from '../modules/filter_clustered_variants.nf'
include { filter_germline_variants }  from '../modules/filter_germline_variants.nf'

workflow preprocess_and_filter_poe {


    take:
        merged_input
        panel_of_errors_fa
        ref
        segdup_regions
        centromere_regions
        simple_repeat_regions
        kg_indels
        regions_input
        germline_inputs

    main: 

    filter_germline_variants(merged_input.join(germline_inputs),regions_input)
    filter_clustered_variants(filter_germline_variants.out.vcf,regions_input)
    filter_centromere_segdups(filter_clustered_variants.out.vcf,segdup_regions,centromere_regions,simple_repeat_regions,kg_indels,regions_input)
    filter_poe(filter_centromere_segdups.out.vcf,panel_of_errors_fa,regions_input)

    emit:
    filter_poe.out.vcf
}
