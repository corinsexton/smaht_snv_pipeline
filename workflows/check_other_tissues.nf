nextflow.enable.dsl=2

include { split_vcf }                       from '../modules/split_vcf.nf'
include { run_minipileup_sr_only_parallel } from '../modules/run_minipileup_sr_only_parallel.nf'
include { merge_minipileup_sr_only_chunks } from '../modules/merge_minipileup_sr_only_chunks.nf'
include { run_minipileup_sr_only }          from '../modules/run_minipileup_sr_only.nf'

workflow check_other_tissues {

    take:
        vcf_inputs
        ref_input
        bam_inputs
        regions_input
        genome_chunks

    main:

    vcf_inputs
      .join(bam_inputs)
      .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, tissue_ids ->
        tuple(id, vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, tissue_ids)
      }
      .set { vcf_bam_channel }

    // Step 1: Split VCF into genome chunks
    ch_vcf_for_split = vcf_inputs
        .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi ->
               tuple(id, vcf, tbi) }

    ch_split_vcf_out = split_vcf(ch_vcf_for_split, genome_chunks)

    ch_split_vcf_out
        .flatMap { id, vcfs, tbis ->
            assert vcfs.size() == tbis.size(), "VCF/TBI list length mismatch for ${id}"
            vcfs.indices.collect { i -> tuple(id, vcfs[i], tbis[i]) }
        }
        .set { ch_split_vcf_chunks }

    // Step 2: Rejoin split VCF chunks with BAM metadata
    ch_split_vcf_chunks
        .combine(vcf_bam_channel, by: 0)
        .map { id, chunk_vcf, chunk_tbi, full_vcf, full_tbi, truth_vcf, truth_vcf_tbi,
                sr_bams, sr_bais, sr_ids ->
                tuple(id, chunk_vcf, chunk_tbi, truth_vcf, truth_vcf_tbi,
                      sr_bams, sr_bais, sr_ids) }
        .set { vcf_chunk_metadata }

    // Step 3: Parallel SR-only minipileup on each chunk
    run_minipileup_sr_only_parallel(vcf_chunk_metadata, ref_input, regions_input)

    // Step 4: Group all chunk outputs per id
    run_minipileup_sr_only_parallel.out.vcf
        .groupTuple(size: 52)
        .map { id, chunk_vcfs, chunk_tbis, truth_vcfs, truth_tbis, mp_vcfs, mp_tbis ->
                tuple(id, mp_vcfs, mp_tbis, truth_vcfs.unique(), truth_tbis.unique()) }
        .set { chunk_groups }

    // Step 5: Join chunk groups with original full VCF and merge
    chunk_groups
        .join(vcf_bam_channel)
        .map { id, mp_vcfs, mp_tbis, truth_vcfs, truth_tbis,
                vcf, tbi, truth_vcf, truth_vcf_tbi, sr_bams, sr_bais, sr_ids ->
                tuple(id, mp_vcfs, mp_tbis, vcf, tbi, truth_vcf, truth_vcf_tbi) }
        .set { merged_sr_input }

    merge_minipileup_sr_only_chunks(merged_sr_input)

    // Step 6: Post-process (parse cross-tissue, set filter, metrics)
    run_minipileup_sr_only(merge_minipileup_sr_only_chunks.out, ref_input, regions_input)

    emit:
    run_minipileup_sr_only.out.vcf
}
