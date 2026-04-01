# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Pipeline Does

This is a Nextflow DSL2 pipeline for calling and filtering mosaic somatic SNVs in the SMaHT consortium. It merges variants from multiple callers (TNhaplotyper2, Strelka2, longcallD, RUFUS), applies sequential filters, annotates with VEP, tiers variants by long-read support, and validates phasing and cross-tissue presence.

## Running the Pipeline

The pipeline is submitted to SLURM via a wrapper script:

```bash
# Edit run_snv_pipeline.sh with the appropriate CSV paths and results dir, then:
sbatch run_snv_pipeline.sh
```

Direct Nextflow invocation example:
```bash
nextflow run main.nf -resume \
  --vep_config vep.ini \
  --longread_csv <lr.csv> \
  --ont_csv <ont.csv> \
  --shortread_csv <sr.csv> \
  --input_metadata <metadata.csv> \
  --input_vcfs <vcfs.csv> \
  --results_dir ./results
```

Key flags: `-resume` restores cached intermediate results; `-with-report` generates an HTML execution report.

## Input CSV Formats

| CSV param | Required columns |
|-----------|-----------------|
| `--input_vcfs` | `id,caller,vcf,vcf.tbi` |
| `--shortread_csv` | `id,cram,crai` |
| `--longread_csv` | `id,cram,crai` |
| `--ont_csv` | `id,cram,crai` (optional) |
| `--input_metadata` | `id,truth_vcf,germline_vcf,sex` |

IDs must match across CSVs. Multiple rows per donor/tissue are grouped internally.

## Architecture Overview

The pipeline is structured into three layers:

**`main.nf`** — Entry point. Parses input CSVs, groups data by donor/tissue, and chains the 6 sub-workflows sequentially.

**`workflows/`** — Six sub-workflows that compose module processes into logical stages:

| Workflow | Stage | What it does |
|----------|-------|-------------|
| `preprocess_merge_callers.nf` | 1–2 | Normalize/PASS-filter VCFs; merge 3–4 callers |
| `preprocess_and_filter_poe.nf` | 3–6 | Remove germline, clustered, PoE, segdup/centromere variants |
| `run_vep.nf` | 7 | Annotate with VEP + gnomAD 4.1 AFs |
| `split_tier1_tier2.nf` | 8–10 | Minipileup SR+LR, binomial tiering (Tier1=LR-supported, Tier2=LR-absent) |
| `phasing.nf` (via `preprocess_and_filter_poe.nf`) | 11 | Phase mosaic variants using germline haplotypes + LR reads |
| `check_other_tissues.nf` | 12–13 | Cross-tissue SR minipileup and final filtering |

**`modules/`** — Individual Nextflow processes (one task per file). Each calls tools or scripts in `bin/`.

**`bin/`** — Python and Bash helper scripts implementing the actual algorithms:
- `merge_callers.py` — Multi-caller VCF reconciliation
- `filter_by_poe.py` — IUPAC-based Panel of Errors matching
- `filter_clustered_variants.py` — Remove variants within proximity threshold
- `tier_filter_variants_SR_PB_ONT.py` — Binomial test tiering using SR/LR/ONT read counts
- `phasing_step2_phase_mosaic.py` — Binomial phasing against germline haplotypes
- `set_filter_vcf.py` — Set final FILTER field
- `bcftools_PASS_norm_dedup.sh` — Normalize, PASS-filter, deduplicate

## Data Flow

```
Input VCFs (4 callers per sample)
  → preprocess_merge_callers   → merged VCF per sample
  → preprocess_and_filter_poe  → filtered VCF (germline/PoE/clustered/regions removed)
  → run_vep                    → VEP-annotated VCF
  → split_tier1_tier2          → tiered VCF (HighConf/LowConf tags)
  → phasing                    → phased VCF with phase tags
  → check_other_tissues        → final VCF with cross-tissue support flags
```

Each stage writes output to a numbered subdirectory of `--results_dir` (e.g., `1_pass_filtered/`, `2_merged_vcf/`, ..., `13_final/`) along with TSV metrics files tracking variant counts before/after filtering.

## Execution Environment

- **Executor**: SLURM (`nextflow.config`)
- **Containers**: Singularity; bcftools uses a custom container built by cos689; VEP uses `park-smaht_dac` container with gnomAD 4.1 custom annotation
- **Resources**: Dynamic queue selection (short/medium/park) based on CPU/memory/time; max 64 GB RAM, 36 CPUs, 48h walltime, 3 retries on error
- **Reference genome**: GRCh38 (no alt contigs), paths defined as params in `main.nf`

## Key Reference Data Paths (in `main.nf` params)

- Reference FASTA: GRCh38 no-alt
- Panel of Errors: `PON.q20q20.05.5.fa.gz`
- Region BEDs: easy/difficult/extreme SMaHT region definitions
- Segdup/centromere BED
- Mills & 1000G indel set (for BQSR context)
