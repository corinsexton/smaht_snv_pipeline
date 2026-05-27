#!/usr/bin/env nextflow


nextflow.enable.dsl=2

include { preprocess_merge_callers } from './workflows/merge_callers.nf'
include { preprocess_and_filter_poe } from './workflows/preprocess_and_filter_poe.nf'
include { run_vep } from './workflows/run_vep'
include { split_tier1_tier2 } from './workflows/split_tier1_tier2.nf'
include { phasing } from './workflows/phasing.nf'
include { check_other_tissues } from './workflows/check_other_tissues.nf'

params.genome_chunks   = "/n/data1/hms/dbmi/park/corinne/smaht/smahtSNV_v2_core_specific/smaht_snv_pipeline/conf/genome_chunks_chr.txt"
params.longread_csv    = null
params.ont_csv         = null
params.panel_of_errors ="/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/smaht_snv_pipeline/panel_of_errors/PON.q20q20.05.5.fa.gz"
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

def poe_input = file(params.panel_of_errors)

def regions_input = tuple(file(params.easy_regions),file(params.diff_regions),file(params.ext_regions),
                          file(params.easy_regions + ".tbi"),file(params.diff_regions + ".tbi"),file(params.ext_regions + ".tbi"))

def ref_fa = file(params.ref)
def ref_fai = file(params.ref_index)
def ref_dict = file(params.ref_dict)
def ref_input = tuple(ref_fa, ref_fai, ref_dict)


// ---------- helper to parse CRAM samplesheet ----------
// SR/LR format: tissue,core,cram,crai (header-based, one row per CRAM file)
def parse_cram_csv(csv_path) {
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.tissue.trim(),
                row.core.trim(),
                file(row.cram.trim()),
                file(row.crai.trim())
            )
        }
}

// ONT format: tissue,cram,crai (no core — ONT is always used pooled per donor)
def parse_pooled_cram_csv(csv_path) {
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.tissue.trim(),
                file(row.cram.trim()),
                file(row.crai.trim())
            )
        }
}

def input_sr  = parse_cram_csv(params.shortread_csv)
def input_lr  = params.longread_csv ? parse_cram_csv(params.longread_csv)        : Channel.empty()
def input_ont = params.ont_csv      ? parse_pooled_cram_csv(params.ont_csv)      : Channel.empty()

// ---------- SR: tissue-level pool (all cores pooled, for minipileup) ----------
input_sr
    .map { tissue, core, cram, crai -> tuple(tissue, cram, crai) }
    .groupTuple(by: 0)
    .map { tissue, crams, crais -> tuple(tissue, crams, crais) }
    .set { sr_by_tissue }

// ---------- SR: donor-level pool (all tissues, for cross-tissue check) ----------
input_sr
    .map { tissue, core, cram, crai ->
        def donor = tissue.tokenize('-')[0]
        tuple(donor, tissue, cram, crai)
    }
    .groupTuple(by: 0)
    .map { donor, tissues, crams, crais ->
        tuple(donor, crams, crais, tissues)   // tissue label per CRAM = source tissue
    }
    .flatMap { donor, crams, crais, tissue_labels ->
        tissue_labels.unique(false).collect { tissue ->
            tuple(tissue, crams, crais, tissue_labels)
        }
    }
    .set { sr_by_donor }

// ---------- helper: donor → tissue mapping (for LR/ONT projection) ----------
sr_by_tissue
    .map { tissue, crams, crais ->
        def donor = tissue.tokenize('-')[0]
        tuple(donor, tissue)
    }
    .set { sr_ids_by_donor }

// ---------- LR: donor-level aggregate (optional; map lookup handles missing donors) ----------
input_lr
    .map { tissue, core, cram, crai ->
        def donor = tissue.tokenize('-')[0]
        tuple(donor, tissue, cram, crai)
    }
    .groupTuple(by: 0)
    .map { donor, tissues, crams, crais ->
        tuple(donor, crams, crais, tissues)
    }
    .set { lr_donor_agg }

// Project LR onto every SR tissue for that donor (map lookup handles missing donors)
lr_donor_map = lr_donor_agg
    .map { donor, crams, crais, tissues ->
        tuple(donor, tuple(crams, crais, tissues))
    }
    .toList()
    .map { list -> list.collectEntries() }

sr_ids_by_donor
    .combine(lr_donor_map)
    .map { donor, tissue, lr_map ->
        def lr_data    = lr_map[donor]
        def lr_crams   = lr_data ? lr_data[0] : []
        def lr_crais   = lr_data ? lr_data[1] : []
        def lr_tissues = lr_data ? lr_data[2] : []
        tuple(tissue, lr_crams, lr_crais, lr_tissues)
    }
    .set { lr_by_tissue }

// ---------- ONT: donor-level aggregate (optional; map lookup handles missing donors) ----------
input_ont
    .map { tissue, cram, crai ->
        def donor = tissue.tokenize('-')[0]
        tuple(donor, tissue, cram, crai)
    }
    .groupTuple(by: 0)
    .map { donor, tissues, crams, crais ->
        tuple(donor, crams, crais, tissues)
    }
    .set { ont_donor_agg }

ont_donor_map = ont_donor_agg
    .map { donor, crams, crais, tissues ->
        tuple(donor, tuple(crams, crais, tissues))
    }
    .toList()
    .map { list -> list.collectEntries() }

sr_ids_by_donor
    .combine(ont_donor_map)
    .map { donor, tissue, ont_map ->
        def ont_data    = ont_map[donor]
        def ont_crams   = ont_data ? ont_data[0] : []
        def ont_crais   = ont_data ? ont_data[1] : []
        def ont_tissues = ont_data ? ont_data[2] : []
        tuple(tissue, ont_crams, ont_crais, ont_tissues)
    }
    .set { ont_by_tissue }

// ---------- Combined tissue-level BAM channel ----------
def input_bams = sr_by_tissue
    .join(lr_by_tissue)
    .join(ont_by_tissue)
    .map { tissue, sr_crams, sr_crais, lr_crams, lr_crais, lr_tissues, ont_crams, ont_crais, ont_tissues ->
        if (!sr_crams || sr_crams.size() == 0) {
            throw new IllegalArgumentException("Tissue ${tissue} has no short-read CRAMs — at least one required")
        }
        tuple(tissue, sr_crams, sr_crais, lr_crams, lr_crais, lr_tissues, ont_crams, ont_crais, ont_tissues)
    }

// ---------- Core → CRAM-basename mapping per tissue (for tier script) ----------
// Format: core \t cram_basename \t type  (type = SR | PB | ONT)
// Only tissue-matched CRAMs are included (LR/ONT keyed by their source tissue).
def cram_map_sr = input_sr
    .map { tissue, core, cram, crai ->
        def basename = cram.name.replaceAll(/\.(cram|bam)$/, '')
        ["${tissue}.core_cram_map.tsv", "${core}\t${basename}\tSR\n"]
    }

def cram_map_lr = params.longread_csv ? input_lr
    .map { tissue, core, cram, crai ->
        def basename = cram.name.replaceAll(/\.(cram|bam)$/, '')
        ["${tissue}.core_cram_map.tsv", "${core}\t${basename}\tPB\n"]
    } : Channel.empty()

// ONT has no core — not added to core_cram_map (pooled annotation only)

// ---------- Auto-expand merged cores (e.g., 001C1-001A3) ----------
// Merged cores are identified by a dash in the core name in the VCF samplesheet.
// Their constituent individual core CRAMs are found automatically, so merged CRAM
// files do not need to be listed in the CRAM samplesheets.
//
// For each merged core, we look up each constituent individual core in the CRAM
// samplesheet (joining on [tissue, indiv_core]) and add the CRAM basename to the
// core_cram_map under the merged core name.

// Helper closure: parse VCF samplesheet and emit ([tissue, indiv_core], merged_core)
// for every constituent of every merged core found.
def mergedCoreExpansions = {
    Channel
        .fromPath(params.input_vcfs)
        .splitCsv(header: true)
        .map { row -> tuple(row.tissue.trim(), row.core.trim()) }
        .filter { tissue, core -> core.contains('-') }
        .unique()
        .flatMap { tissue, merged_core ->
            merged_core.tokenize('-').collect { indiv_core ->
                tuple([tissue, indiv_core], merged_core)
            }
        }
}

def cram_map_sr_merged = mergedCoreExpansions()
    .join(input_sr.map { tissue, core, cram, crai -> tuple([tissue, core], cram) })
    .map { key, merged_core, cram ->
        def basename = cram.name.replaceAll(/\.(cram|bam)$/, '')
        ["${key[0]}.core_cram_map.tsv", "${merged_core}\t${basename}\tSR\n"]
    }

def cram_map_lr_merged = params.longread_csv ? mergedCoreExpansions()
    .join(input_lr.map { tissue, core, cram, crai -> tuple([tissue, core], cram) })
    .map { key, merged_core, cram ->
        def basename = cram.name.replaceAll(/\.(cram|bam)$/, '')
        ["${key[0]}.core_cram_map.tsv", "${merged_core}\t${basename}\tPB\n"]
    } : Channel.empty()

cram_map_sr.mix(cram_map_lr)
    .mix(cram_map_sr_merged).mix(cram_map_lr_merged)
    .collectFile { item -> item }
    .map { f ->
        def tissue = f.name.replace('.core_cram_map.tsv', '')
        tuple(tissue, f)
    }
    .set { core_cram_map }

/////

// ---------- VCF samplesheet ----------
// New format: tissue,core,gcc,caller,caller_vcf (header-based)
def input_vcfs = Channel
    .fromPath(params.input_vcfs)
    .splitCsv(header: true)
    .map { row ->
        tuple(
            row.tissue.trim(),
            row.core.trim(),
            row.gcc.trim(),
            row.caller.trim(),
            file(row.caller_vcf.trim()),
            file(row.caller_vcf.trim() + '.tbi')
        )
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
        def tbi_path = file(row.germline_calls + '.tbi')
        def germline_tbi = tbi_path.exists() ? tbi_path : file(row.germline_calls + '.csi')
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

    merged_calls = preprocess_merge_callers(
        input_vcfs.combine(truth_ch, by: 0),
        params.ref,
        params.ref_index,
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

    vep_snvs_out = run_vep(vep_input, regions_input)

    tier_split_output = split_tier1_tier2(vep_snvs_out.join(truth_ch), input_bams, ref_input, regions_input, file(params.genome_chunks), core_cram_map)

    def has_lr = params.longread_csv || params.ont_csv

    def pipeline_output = has_lr
        ? phasing(tier_split_output, germline_calls_ch, input_bams, ref_input, vep_config, regions_input, sex_ch)
        : tier_split_output

    check_other_tissues(pipeline_output, ref_input, sr_by_donor, regions_input, file(params.genome_chunks))

}

