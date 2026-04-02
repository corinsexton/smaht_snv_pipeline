/*
 * Workflow to preprocess VCFs and merge
 *
 * Steps:
 *  1. Normalize, filter PASS, and atomize VCF (preprocess_vcf)
 *  2. Collect all processed VCFs per tissue (all cores, all callers)
 *  3. Merge callers for each tissue
 */

nextflow.enable.dsl=2

include { preprocess_vcf }  from '../modules/preprocess_vcf.nf'
include { merge_callers }   from '../modules/merge_callers.nf'

workflow preprocess_merge_callers {

    take:
        vcf_inputs      // (tissue, core, gcc, caller, vcf, tbi, truth, truth_tbi) — one row per caller VCF
        ref
        ref_index
        regions_input

    main:
        // Step 1: Normalize, filter PASS, atomize
        // Core is embedded in the caller field as CORE__CALLER so it survives preprocess_vcf
        // (preprocess_vcf uses caller in output filename; __ is safe for filesystems)
        preprocess_input = vcf_inputs
            .map { tissue, core, gcc, caller, vcf, tbi, truth, truth_tbi ->
                tuple(tissue, "${core}__${caller}", vcf, tbi, ref, ref_index, truth, truth_tbi)
            }
        preprocess_vcf(preprocess_input, regions_input)

        // Step 2: Group processed VCFs by tissue
        // core_callers list retains CORE__CALLER encoding for merge step
        grouped = preprocess_vcf.out.vcf
            .groupTuple(by: 0)
            .map { tissue, core_callers, vcfs, tbis, truth_vcfs, truth_tbis ->
                tuple(
                    tissue,
                    core_callers.unique(),
                    vcfs.unique(),
                    tbis.unique(),
                    truth_vcfs.unique(),
                    truth_tbis.unique()
                )
            }

        // Step 3: Merge all callers for the tissue
        merge_callers(grouped)

    emit:
        merge_callers.out.vcf
}
