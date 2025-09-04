#!/usr/bin/env nextflow


nextflow.enable.dsl=2

include { preprocess_and_filter_poe } from './workflows/preprocess_and_filter_poe.nf'
include { run_vep } from './workflows/run_vep'
include { split_tier1_tier2 } from './workflows/split_tier1_tier2.nf'

//params.input_csv       = "./mt_calls.csv"
//params.input_csv       = "./150x_files_ALL.csv"
//params.input_csv       = "./150x_files.csv"
//params.input_csv       = "./150x_files_w_bams.csv"
params.panel_of_errors = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/panel_of_errors/PON.q20q20.05.5.fa.gz"
params.panel_of_errors_index  = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/panel_of_errors/PON.q20q20.05.5.fa.gz.fai"
params.results_dir     = "./new_results"
params.ref             = "/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa"
params.ref_index             = "/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa.fai"
params.segdup_regions = "/n/data1/hms/dbmi/park/SOFTWARE/UCSC/GRCh38_UCSC_SegDup.formatted.bed"
params.centromere_regions = "/home/cos689/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/centromere.bed"
params.max_depth = 300

def ensureTabixIndex(vcf_path) {
    def tbi_path = file("${vcf_path}.tbi")
    //println "Checking VCF: ${vcf_path}"

    if (!vcf_path.exists()) {
        error "VCF not found: ${vcf_path}"
    }

    if (!tbi_path.exists()) {
        //println "Index missing — running tabix: ${vcf_path}"
        def proc = ["tabix", "-p", "vcf", vcf_path.toString()].execute()
        proc.waitFor()
        //println "Tabix exit code: ${proc.exitValue()}"
        if (!tbi_path.exists()) {
            error "Failed to create index: ${tbi_path}"
        }
    }

    return tbi_path
}

def poe_fa = file(params.panel_of_errors)
def poe_fai = file(params.panel_of_errors_index)
def poe_input = tuple(poe_fa, poe_fai)

def ref_fa = file(params.ref)
def ref_fai = file(params.ref_index)
def ref_input = tuple(ref_fa, ref_fai)

def input_vcfs = Channel
    .fromPath(params.input_csv)
    .splitCsv(header: true)
    .map { row ->
        def id  = row.id
        def vcf = file(row.vcf)
        def tbi = ensureTabixIndex(vcf)
        tuple(id, vcf, tbi)
    }

def input_bams = Channel
    .fromPath(params.input_csv)
    .splitCsv(header: true)
    .map { row ->
        def id = row.id
        def bam = file(row.bam)
        def bai = file(row.bai)
        def lr_bam = file(row.lr_bam)
        def lr_bai = file(row.lr_bai)
        tuple(id,bam,bai,lr_bam,lr_bai)
    }

workflow {

    // remove first filters
    filtered = preprocess_and_filter_poe(
        input_vcfs,
        poe_input,
        params.ref,
        params.segdup_regions,
        params.centromere_regions,
        params.max_depth
    )

    // use VEP for AF filtering
    vep_config = Channel.fromPath(params.vep_config)
    filtered.combine( vep_config )
                        .map {
                              id, vcf, tbi, config ->
                               [id: id, file: vcf, index: tbi, vep_config: config]
                             }
                        .set {vep_input}

            
    vep_snvs_out = run_vep(vep_input) // output: id, snv_vcf, snv_tbi

    // run pileup and split based on LR presence (tier1) / absence (tier2)
    tier_split_output = split_tier1_tier2(vep_snvs_out, input_bams, ref_input)


    // run last filters based on tier1 or tier2
    // tier1_filters(tier_split_output.tier1)
    // tier2_filters(tier_split_output.tier2)
    

}

