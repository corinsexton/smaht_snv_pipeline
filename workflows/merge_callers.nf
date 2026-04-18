/*
 * Workflow to preprocess VCFs and merge
 *
 * Steps:
 *  1. Normalize, filter PASS, and atomize VCF (preprocess_vcf)
 *  2. Collect all processed VCFs per ID (3–4 callers)
 *  3. Merge callers for each ID
 */

nextflow.enable.dsl=2

include { preprocess_vcf }       from '../modules/preprocess_vcf.nf'
include { merge_callers }  from '../modules/merge_callers.nf'

workflow preprocess_merge_callers {


    take:
        vcf_inputs
        ref
        ref_index
        regions_input

    main: 
        // Step 1: Normalize, filter PASS, atomize
        // Assumes: preprocess_vcf emits (id, processed_vcf, processed_vcf.tbi)
        preprocess_input = vcf_inputs
            .map { id, core, gcc, caller, vcf, tbi, truth, truth_tbi ->
                tuple(id, core, caller, vcf, tbi, ref, ref_index, truth, truth_tbi)
            }
        preprocess_vcf(preprocess_input, regions_input)

        // ---- Step 2: group processed VCFs by sample ID ----
        grouped = preprocess_vcf.out.vcf
                   .groupTuple(by: 0)
                   .map { id, cores, callers, vcfs, tbis, truth_vcfs, truth_tbis ->
                    def entries = [cores, callers, vcfs, tbis].transpose()
                                    .unique { a -> "${a[0]}__${a[1]}" }
                                    .sort   { a, b -> "${a[0]}__${a[1]}" <=> "${b[0]}__${b[1]}" }
                    tuple(
                        id,
                        entries.collect { c, cal, v, t -> "${c}__${cal}" },
                        entries.collect { c, cal, v, t -> v },
                        entries.collect { c, cal, v, t -> t },
                        truth_vcfs.unique(),
                        truth_tbis.unique()
                    )
                }
                   // grouped emits:
                   // tuple(id, [ [caller1 tuple], [caller2 tuple], ... ])

        // ---- Step 3: send grouped lists to merge_callers ----
        merge_callers(grouped)

    emit:
        merge_callers.out.vcf
}
