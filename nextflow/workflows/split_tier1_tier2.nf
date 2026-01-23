nextflow.enable.dsl=2


include { run_minipileup_parallel } from '../modules/run_minipileup_parallel.nf'
include { tier_variants } from '../modules/tier_variants.nf'
include { filter_binom_fisher } from '../modules/filter_binom_fisher.nf'
include { split_vcf } from '../modules/split_vcf.nf'
include { merge_minipileup_chunks } from '../modules/merge_minipileup_chunks.nf' 

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
      .set { vcf_bam_joined }



     /* * Step 1: Split VCF into 3 or 4 chunks per sample */ 
    ch_vcf_for_split = vcf_inputs 
        .map { id, vcf, tbi, truth_vcf, truth_vcf_tbi -> 
               tuple(id, vcf, tbi) } 
    
    ch_split_vcf_out = split_vcf(ch_vcf_for_split)


    ch_split_vcf_out
        .flatMap { id, vcfs, tbis ->
            assert vcfs.size() == tbis.size(), "VCF/TBI list length mismatch for ${id}"
            vcfs.indices.collect { i -> tuple(id, vcfs[i], tbis[i]) }
        }
        .set { ch_split_vcf_chunks }

    /* * Step 2: Rejoin split VCF chunks with BAM metadata */ 
    ch_split_vcf_chunks 
        .combine(vcf_bam_joined, by: 0) 
        .map { id, chunk_vcf, chunk_tbi, full_vcf, full_tbi, truth_vcf, truth_vcf_tbi,
                sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais -> 
                tuple(id, chunk_vcf, chunk_tbi, truth_vcf, truth_vcf_tbi, 
                sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais) } 
        .set { vcf_chunk_metadata } 

    // Step 3: Parallel minipileup on each chunk 
    run_minipileup_parallel(vcf_chunk_metadata, ref_input) 
    
    // Step 4: Group all chunk outputs per id 
    run_minipileup_parallel.out.vcf 
        .groupTuple(size:24) 
        .map { id, chunk_vcfs, chunk_tbis, truth_vcfs, truth_tbis, mp_vcfs, mp_tbis -> 
                tuple( id, mp_vcfs, mp_tbis, truth_vcfs.unique(), truth_tbis.unique() ) } 
        .set { chunk_groups } 
    
    // Step 5: Join chunk groups with original metadata and merge 
    chunk_groups 
        .join(vcf_bam_joined) 
        .map { id, mp_vcfs, mp_tbis, truth_vcfs, truth_tbis,
                vcf, tbi, truth_vcf, truth_vcf_tbi,
                sr_bams, sr_bais, lr_bams, lr_bais, lr_ont_bams, lr_ont_bais -> 
                tuple(id, mp_vcfs, mp_tbis, vcf, tbi, truth_vcf, truth_vcf_tbi) } 
        .set { merged_minipileups_input } 
    
    merge_minipileup_chunks(merged_minipileups_input) 
    
    // Step 6: Provide single merged minipileup VCF to tier_variants 
    tier_variants( merge_minipileup_chunks.out, regions_input ) 
    
    // Step 7: Downstream filters 
    filter_binom_fisher(tier_variants.out.vcf, regions_input) 
    
    emit: 
        filter_binom_fisher.out.vcf 
}
