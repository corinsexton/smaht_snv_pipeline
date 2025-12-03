/*
 * Workflow to preprocess VCFs and filter using panel of errors (POE)
 *  Filter out known artifact variants from POE VCF (filter_panel_errors)
 */

nextflow.enable.dsl=2

include { preprocess_vcf }       from '../modules/preprocess_vcf.nf'
include { filter_centromere_segdups }  from '../modules/filter_centromere_segdups.nf'
include { filter_poe }  from '../modules/filter_panel_errors.nf'
include { filter_clustered_variants }  from '../modules/filter_clustered_variants.nf'
include { filter_high_cov }  from '../modules/filter_high_cov.nf'

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

    main: 

    filter_clustered_variants(merged_input,regions_input)
    filter_centromere_segdups(filter_clustered_variants.out.vcf,segdup_regions,centromere_regions,simple_repeat_regions,kg_indels,regions_input)
    filter_poe(filter_centromere_segdups.out.vcf,panel_of_errors_fa,regions_input)

    emit:
    filter_poe.out.vcf
}
