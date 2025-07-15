#!/usr/bin/env nextflow


nextflow.enable.dsl=2

include { preprocess_and_filter_poe } from './workflows/preprocess_and_filter_poe.nf'
include { run_vep } from './workflows/run_vep'

//params.input_csv       = "./mt_calls.csv"
params.input_csv       = "./150x_files_ALL.csv"
params.panel_of_errors = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/panel_of_errors/POE_benchmarking.vcf.gz"
params.results_dir     = "./results"
params.ref             = "/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa"
params.segdup_regions = "/n/data1/hms/dbmi/park/SOFTWARE/UCSC/GRCh38_UCSC_SegDup.formatted.bed"
params.centromere_regions = "/home/cos689/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/centromere.bed"

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

def poe_vcf = file(params.panel_of_errors)
def poe_tbi = ensureTabixIndex(poe_vcf)
def poe_input = Channel.of(tuple(poe_vcf, poe_tbi))

def input_vcfs = Channel
    .fromPath(params.input_csv)
    .splitCsv(header: true)
    .map { row ->
        def id  = row.id
        def vcf = file(row.vcf)
        def tbi = ensureTabixIndex(vcf)
        tuple(id, vcf, tbi)
    }

workflow {

    // Run full workflow
    filtered = preprocess_and_filter_poe(
        input_vcfs,
        poe_input,
        params.ref,
        params.segdup_regions,
        params.centromere_regions
    )

    vep_config = Channel.fromPath(params.vep_config)
    filtered.combine( vep_config )
                        .map {
                              id, vcf, tbi, config ->
                               [id: id, file: vcf, index: tbi, vep_config: config]
                             }
                        .set {vep_input}

            
    run_vep(vep_input)

    //// Save final outputs to results/filtered/
    //run_vep
    //    .subscribe { id, vcf, tbi ->
    //        file("${params.results_dir}").mkdirs()
    //        vcf.copyTo("${params.results_dir}")
    //        tbi.copyTo("${params.results_dir}")
    //    }
}

