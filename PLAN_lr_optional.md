# Plan: Make LR Completely Optional

## Branch

All edits below should be made on a new branch named **`no-long-read`**:
```bash
git checkout -b no-long-read
```

---

## Context

The pipeline currently fails when `--longread_csv` and `--ont_csv` are omitted. Two blockers:

1. `minipileup2-parallel.sh` (the external tool called by `run_minipileup2_parallel`) has a hard-coded check: `(( ${#PB_CRAMS[@]} > 0 )) || { echo "Error: at least one PB long-read required"; exit 1; }` — it exits non-zero if no PB CRAMs are passed.
2. `phasing` workflow always runs; `phase_mosaic_vars.sh` requires `--pb-cram`. User confirmed: **skip phasing entirely when no LR is provided**.

The goal: pass only `--shortread_csv` and the pipeline completes SR-only. Providing both LR CSVs preserves the existing LR-augmented path exactly.

---

## What already works (no change needed)

- **Channel construction in `main.nf`** (`lr_donor_map` / `ont_donor_map`): when `input_lr = Channel.empty()`, `toList()` emits `[]`, `collectEntries()` returns `[:]`, and the map-lookup in `sr_ids_by_donor.combine(lr_donor_map)` returns `null` for every donor → propagates `(tissue, [], [], [])` to `lr_by_tissue` / `ont_by_tissue`. The subsequent `.join()` to build `input_bams` produces 9-element tuples with empty lists for LR/ONT fields.
- **`cram_map_lr` / `cram_map_lr_merged`**: already guarded by `params.longread_csv ? ... : Channel.empty()`.
- **Bash for-loops in process scripts**: `for f in ${lr_bams}; do` is a no-op when the staged path list is empty; `pb_lr_args` / `ont_lr_args` remain `""`.

---

## Changes required

### 1. `main.nf` — two edits

**a. Add param defaults** (near the top with the other `params.*` declarations, ~line 13):
```groovy
params.longread_csv = null
params.ont_csv      = null
```
Prevents Nextflow from complaining if a user omits these; makes the intent explicit.

**b. Derive `has_lr` and gate phasing** (in the `workflow {}` block, ~lines 345–349):

Current:
```groovy
tier_split_output = split_tier1_tier2(...)
phasing_output = phasing(tier_split_output, germline_calls_ch, input_bams, ...)
check_other_tissues(phasing_output, ...)
```

Replace with:
```groovy
def has_lr = params.longread_csv || params.ont_csv

tier_split_output = split_tier1_tier2(...)

def pipeline_output = has_lr
    ? phasing(tier_split_output, germline_calls_ch, input_bams, ref_input, vep_config, regions_input, sex_ch)
    : tier_split_output

check_other_tissues(pipeline_output, ref_input, sr_by_donor, regions_input, file(params.genome_chunks))
```

`tier_split_output` and `phasing_output` share the same 5-tuple shape `(id, vcf.gz, vcf.gz.tbi, truth_vcf, truth_vcf_tbi)`, so `check_other_tissues` accepts either without modification.

---

### 2. `modules/run_minipileup2_parallel.nf` — bash script conditional

After the array-building loops (after `ont_lr_args` is set), replace the single `minipileup2-parallel.sh` call with a conditional:

```bash
if [[ ${#pb_cram_arr[@]} -eq 0 && ${#ont_cram_arr[@]} -eq 0 ]]; then
    # No LR available — fall back to SR-only minipileup
    minipileup2-parallel_sr_only.sh -i ${vcf} \
        -r ${ref} \
        -o ${id}.${chr}.minipileup \
        ${sr_crams}
else
    minipileup2-parallel.sh -i ${vcf} \
        -r ${ref} \
        -o ${id}.${chr}.minipileup \
        ${sr_crams} \
        ${pb_lr_args} \
        ${ont_lr_args}
fi
```

**Assumption**: `minipileup2-parallel_sr_only.sh` treats `-o <prefix>` literally (appends `.vcf.gz`, not `_sr.vcf.gz`). This is consistent with how the cross-tissue module uses it — it passes `-o ${id}.${chunk}.minipileup_sr` and the `_sr` is in the prefix string, not added by the tool. The output glob `path("${id}*.minipileup.vcf.gz")` will match either way.

**If that assumption is wrong** (tool adds `_sr` internally): add a rename step after the SR-only call:
```bash
mv ${id}.${chr}.minipileup_sr.vcf.gz    ${id}.${chr}.minipileup.vcf.gz
mv ${id}.${chr}.minipileup_sr.vcf.gz.tbi ${id}.${chr}.minipileup.vcf.gz.tbi
```

---

### 3. `run_snv_pipeline.sh` — remove hardcoded LR params

Remove (or comment out) the two LR lines so the script runs SR-only by default:
```bash
# --longread_csv samplesheets/p25_lr.csv \
# --ont_csv samplesheets/p25_ont.csv \
```
When you do want to run with LR, add them back (or keep a separate run script for the LR case).

---

## Files NOT changing

| File | Reason |
|------|--------|
| `workflows/split_tier1_tier2.nf` | LR fields already flow through as empty lists; no branching needed at this level |
| `workflows/phasing.nf` | Now only called when `has_lr` is true; no internal changes |
| `workflows/check_other_tissues.nf` | Already SR-only; takes 5-tuple input compatible with both `tier_split_output` and `phasing_output` |
| `nextflow.config` | No resource or param changes needed |
| `main.nf` channel construction | Already handles `Channel.empty()` correctly |

---

## Verification

1. **SR-only run**: Execute with only `--shortread_csv` (no LR params). Confirm:
   - `run_minipileup2_parallel` jobs complete without "at least one PB long-read required" error
   - `tier_variants_binom` produces `.tiered.vcf.gz` output
   - Phasing is skipped; `check_other_tissues` receives `tier_split_output` directly
   - Final VCFs land in `results_dir/`

2. **LR run regression**: Execute with all three CSVs. Confirm existing LR path is unchanged — phasing still runs, minipileup2 still uses PB/ONT args.

3. **Mixed-donor run**: If some donors have LR and some don't, the map-lookup already handles per-donor absence; this is unchanged.

4. **PHASING tags in final VCF (SR-only path)**: When phasing is skipped, `check_other_tissues` receives `tier_split_output` directly instead of `phasing_output`. Confirm the final VCF produced in `results_dir/` is well-formed:
   - No dangling `FORMAT` or `INFO` fields that reference phasing tags (e.g. `PS`, `PID`, `PGT`) — these would appear in the header but have no corresponding values if phasing was partially wired. Run `bcftools stats` and check for header/record mismatches.
   - If `tier_split_output` VCFs are produced by a step that conditionally adds phasing `FORMAT` tags when phasing ran, verify those tags are absent (not blank/`.`) in the SR-only output.
   - Spot-check one output VCF: `bcftools view -h <final.vcf.gz> | grep -E "^##FORMAT.*PS|PGT|PID"` should return nothing on the SR-only path.
