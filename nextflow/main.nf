#!/usr/bin/env nextflow


nextflow.enable.dsl=2

include { preprocess_merge_callers } from './workflows/merge_callers.nf'
include { preprocess_and_filter_poe } from './workflows/preprocess_and_filter_poe.nf'
include { run_vep } from './workflows/run_vep'
include { split_tier1_tier2 } from './workflows/split_tier1_tier2.nf'
include { phasing } from './workflows/phasing.nf'
include { check_other_tissues } from './workflows/check_other_tissues.nf'

params.panel_of_errors = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/panel_of_errors/PON.q20q20.05.5.fa.gz"
params.panel_of_errors_index  = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/panel_of_errors/PON.q20q20.05.5.fa.gz.fai"
params.results_dir     = "./new_results"
params.ref             = "/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa"
params.ref_index             = "/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa.fai"
params.ref_dict             = "/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.dict"
params.segdup_regions = "/home/cos689/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/segdup_GRCh38_official.bed.gz"
params.centromere_regions = "/home/cos689/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/centromeres_GRCh38_official.bed.gz"
params.simple_repeats = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/simple_repeats.bed"
params.kg_indels = "/n/data1/hms/dbmi/park/corinne/ref/mills_1kg_gold_standard_indels.vcf.gz"

params.easy_regions = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/SMaHT_easy_v2.bed.gz"
params.diff_regions = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/SMaHT_difficult_v2.bed.gz"
params.ext_regions = "/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/nextflow/SMaHT_extreme_v2.bed.gz"


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

def regions_input = tuple(file(params.easy_regions),file(params.diff_regions),file(params.ext_regions),
                          file(params.easy_regions + ".tbi"),file(params.diff_regions + ".tbi"),file(params.ext_regions + ".tbi"))

def ref_fa = file(params.ref)
def ref_fai = file(params.ref_index)
def ref_dict = file(params.ref_dict)
def ref_input = tuple(ref_fa, ref_fai, ref_dict)


// ---------- helper to parse cram/crai CSV ----------
def parse_cram_csv(csv_path) {
    Channel
        .fromPath(csv_path)
        .splitCsv(header: false)
        .filter { row ->
            // Skip header line if first column starts with 'id' (case-insensitive)
            !(row[0]?.toString()?.toLowerCase()?.startsWith('id'))
        }
        .map { row ->
            // row is now a List of values, not a Map
            def id = row[0].toString().trim()
            def fields = row.drop(1).findAll { it && it.trim() }

            // Sanity check: ensure pairs of CRAM and CRAI
            if (fields.size() % 2 != 0) {
                log.warn "Row for ${id} has an odd number of CRAM/CRAI entries in ${csv_path} (${fields.size()})"
            }

            def crams = (0..<fields.size()/2).collect { i -> file(fields[2*i].trim()) }
            def crais = (0..<fields.size()/2).collect { i -> file(fields[2*i+1].trim()) }

            //log.info "Parsed ${crams.size()} CRAM/CRAI pairs for ${id}"
            tuple(id, crams, crais)
        }
}

def input_sr  = parse_cram_csv(params.shortread_csv)
def input_lr  = params.longread_csv ? parse_cram_csv(params.longread_csv) : Channel.empty()
def input_ont = params.ont_csv      ? parse_cram_csv(params.ont_csv) : Channel.empty()

// ---------- get all donor sr for final cross tissue check --------
//
// Extract donor + tissue, group by donor, and include CRAMs + CRAIs
//
input_sr
    .map { id, crams, crais ->
        def (donor, tissue) = id.tokenize('-')
        tuple(id, donor, tissue, crams, crais)
    }
    .groupTuple(by: 1)   // group by donor
    .map { ids, donor, tissues, crams, crais ->
        // grouped_entries contains tuples: [id, donor, tissue, crams, crais]

        // Flatten donor-level crams & crais (maintains pairing)
        def flat_crams  = crams.flatten()
        def flat_crais  = crais.flatten()

        // Expand tissue labels so each CRAM has one
        def expanded_tissues = []
        tissues.eachWithIndex { tissue, i ->
            expanded_tissues.addAll( Collections.nCopies(crams[i].size(), tissue) )
        }

        // Return donor-level aggregate
        tuple(donor, ids, flat_crams, flat_crais, expanded_tissues)
    }
    .flatMap { donor, all_ids, crams, crais, tissues ->
        // For each original id, emit donor-level data
        all_ids.collect { id ->
            tuple(id, crams, crais, tissues)
        }
    }
    .set { sr_by_donor }

// all SR ids, keyed by donor, so we can “project” LR donor aggregates onto every SR id
input_sr
    .map { id, crams, crais ->
        def (donor, tissue) = id.tokenize('-')
        tuple(donor, id, tissue)   // key=donor, keep id and tissue
    }
    .set { sr_ids_by_donor }
//sr_ids_by_donor.view()
// ---------- get all donor lr for minipileup --------
//
// Extract donor + tissue, group by donor, and include CRAMs + CRAIs
//
input_lr
    .map { id, crams, crais ->
        def (donor, tissue) = id.tokenize('-')
        tuple(id, donor, crams, crais)
    }
    .groupTuple(by: 1)   // group by donor
    .map { ids, donor, crams, crais ->
        // grouped_entries contains tuples: [id, donor, crams, crais]

        // Flatten donor-level crams & crais (maintains pairing)
        def flat_crams  = crams.flatten()
        def flat_crais  = crais.flatten()

        // Expand tissue labels so each CRAM has one
        def expanded_ids = []
        ids.eachWithIndex { id, i ->
            expanded_ids.addAll( Collections.nCopies(crams[i].size(), id) )
        }

        // Return donor-level aggregate
        tuple(donor, flat_crams, flat_crais, expanded_ids)
    }
    .set { lr_donor_agg }

//lr_donor_agg.view()

// Project LR donor aggregate onto every SR id for that donor
sr_ids_by_donor
    .combine(lr_donor_agg, by: 0) 
    .map { donor, id, tissue, lr_crams, lr_crais, lr_tissues ->
        tuple(id, lr_crams, lr_crais, lr_tissues)
    }
    .set { lr_by_donor }

// ---------- get all donor ont for minipileup --------
//
// Extract donor + tissue, group by donor, and include CRAMs + CRAIs
//
input_ont
    .map { id, crams, crais ->
        def (donor, tissue) = id.tokenize('-')
        tuple(id, donor, crams, crais)
    }
    .groupTuple(by: 1)   // group by donor
    .map { ids, donor, crams, crais ->
        // grouped_entries contains tuples: [id, donor, tissue, crams, crais]

        // Flatten donor-level crams & crais (maintains pairing)
        def flat_crams  = crams.flatten()
        def flat_crais  = crais.flatten()

        // Expand tissue labels so each CRAM has one
        def expanded_ids = []
        ids.eachWithIndex { id, i ->
            expanded_ids.addAll( Collections.nCopies(crams[i].size(), id) )
        }

        // Return donor-level aggregate
        tuple(donor, flat_crams, flat_crais, expanded_ids)
    }
    .set { ont_donor_agg }

// Project LR donor aggregate onto every SR id for that donor
sr_ids_by_donor
    .combine(ont_donor_agg, by: 0) 
    .map { donor, id, tissue, lr_crams, lr_crais, lr_tissues ->
        tuple(id, lr_crams, lr_crais, lr_tissues)
    }
    .set { ont_by_donor }

// ---------- merge all by sample ID ----------
def input_bams = input_sr
    .join(lr_by_donor)
    .join(ont_by_donor)
    .map { id, sr_crams, sr_crais, lr_crams, lr_crais, lr_tissues, ont_crams, ont_crais, ont_tissues ->
            if( !sr_crams || sr_crams.size() == 0 ) {
            throw new IllegalArgumentException("Sample ${id} has no short-read CRAMs — at least one required")
            }

        tuple(id, sr_crams, sr_crais, lr_crams, lr_crais, lr_tissues, ont_crams, ont_crais, ont_tissues)
    }

/////

def input_vcfs = Channel
    .fromPath(params.input_vcfs)
    .splitCsv(header: true)
    .map { row ->
        def id  = row.id
        def caller = row.caller
        def vcf = file(row.vcf)
        def tbi = file(row.vcf + '.tbi')
        tuple(id, caller, vcf, tbi)
    }

def truth_ch = Channel
    .fromPath(params.input_metadata)
    .splitCsv(header: true)
    .map{ row ->
        def id = row.id
        def truth_vcf = row.truth_vcf ? file(row.truth_vcf) : file("none.vcf")
        def truth_tbi = row.truth_vcf ? file(row.truth_vcf + '.tbi') : file("none.vcf.tbi")
        tuple(id, truth_vcf, truth_tbi)
    }

def germline_calls_ch = Channel
    .fromPath(params.input_metadata)
    .splitCsv(header: true)
    .map{ row ->
        def id = row.id
        def germline_vcf = row.germline_calls
        def germline_tbi = file(row.germline_calls + '.tbi')
        tuple(id, germline_vcf, germline_tbi)
    }

def sex_ch = Channel
    .fromPath(params.input_metadata)
    .splitCsv(header: true)
    .map{ row ->
        def id = row.id
        def sex = row.sex
        tuple(id, sex)
    }


workflow {

    // remove first filters
    merged_calls = preprocess_merge_callers(
        input_vcfs.combine(truth_ch, by: 0),
        params.ref,
        regions_input
    )

    filtered = preprocess_and_filter_poe(
        merged_calls,
        poe_input,
        params.ref,
        params.segdup_regions,
        params.centromere_regions,
        params.simple_repeats,
        params.kg_indels,
        regions_input,
        germline_calls_ch
    )

    // use VEP for AF filtering
    vep_config = Channel.fromPath(params.vep_config)
    filtered.combine( vep_config )
                        .map {
                              id, vcf, tbi, truth_vcf, truth_tbi, config ->
                               [id: id, file: vcf, index: tbi, truth_vcf:truth_vcf, truth_tbi:truth_tbi, vep_config: config]
                             }
                        .set {vep_input}

            
    vep_snvs_out = run_vep(vep_input, regions_input) // output: id, snv_vcf, snv_tbi

    // run pileup and split based on LR presence (tier1) / absence (tier2)
    tier_split_output = split_tier1_tier2(vep_snvs_out.join(truth_ch), input_bams, ref_input, regions_input)

    phasing_output = phasing(tier_split_output, germline_calls_ch, input_bams, ref_input, vep_config, regions_input, sex_ch)

    check_other_tissues(phasing_output, ref_input, sr_by_donor, regions_input)

}

