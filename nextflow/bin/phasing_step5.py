#!/usr/bin/env python3
import pysam
import csv
import argparse
from collections import Counter

# ---------------------------------------------------------------
# Step 5–7 Combined Script (Sex-aware haplotype classification + phasing tag output)
#   1. Extract reads covering both somatic and germline sites from BAMs
#   2. Count read combinations (both, germ_only, var_only, none)
#   3. Classify haplotypes using defined logic
#   4. Output phasing_tags.tsv for bcftools annotation
# ---------------------------------------------------------------

def get_base_at_pos(aln, chrom, pos):
    """Return base in read at 1-based genome position (None if not covered)."""
    for qpos, refpos in aln.get_aligned_pairs(matches_only=True):
        if refpos == pos - 1:
            return aln.query_sequence[qpos]
    return None

def count_haplotypes(tsv_file, bam_list, sex):
    bam_files = [pysam.AlignmentFile(bam, "rb") for bam in bam_list]
    results = []

    with open(tsv_file) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        for row in reader:
            chrom = row["chrom"]
            var_pos = int(row["Var_pos"])
            germ_pos = row["Germ_pos"]
            if germ_pos == "NA":
                row.update({k: "0" for k in ["case_both","case_germ_only","case_var_only","case_none","total_reads"]})
                row["hap_classification"] = "no_germline"
                results.append(row)
                continue

            germ_pos = int(germ_pos)
            var_alt = row["Var_alt"]
            germ_alt = row["Germ_alt"]

            counts = Counter()
            for bam in bam_files:
                for read in bam.fetch(chrom, min(var_pos, germ_pos)-1, max(var_pos, germ_pos)):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    base_var = get_base_at_pos(read, chrom, var_pos)
                    base_germ = get_base_at_pos(read, chrom, germ_pos)
                    if base_var is None or base_germ is None:
                        continue

                    if base_var == var_alt and base_germ == germ_alt:
                        counts["case_both"] += 1
                    elif base_var == var_alt and base_germ != germ_alt:
                        counts["case_var_only"] += 1
                    elif base_var != var_alt and base_germ == germ_alt:
                        counts["case_germ_only"] += 1
                    else:
                        counts["case_none"] += 1

            total_reads = sum(counts.values())
            for key in ["case_both","case_germ_only","case_var_only","case_none"]:
                row[key] = counts.get(key, 0)
            row["total_reads"] = total_reads

            case_both = counts.get("case_both", 0)
            case_none = counts.get("case_none", 0)
            case_var_only = counts.get("case_var_only", 0)
            case_germ_only = counts.get("case_germ_only", 0)

            # ------------------------------
            # Sex-specific chromosome logic
            # ------------------------------
            if sex == "male":
                # Haploid X/Y in males
                if chrom in ["chrX", "chrY"]:
                    if total_reads == 0:
                        row["hap_classification"] = "no_coverage"
                    else:
                        row["hap_classification"] = "haploid"
                    results.append(row)
                    continue

            elif sex == "female":
                # Females: X diploid, Y = artifact
                if chrom == "chrY":
                    row["hap_classification"] = "artifact"
                    results.append(row)
                    continue
                # chrX treated as diploid below

            else:
                # Unknown sex: treat X/Y as artifacts
                if chrom in ["chrX", "chrY"]:
                    row["hap_classification"] = "artifact"
                    results.append(row)
                    continue

            # ------------------------------
            # Diploid classification logic
            # ------------------------------
            if case_both > 1 and case_none > 1 and case_var_only <= 1 and case_germ_only <= 1:
                row["hap_classification"] = "germline"
            elif case_both <= 1 and case_none <= 1 and case_var_only > 1 and case_germ_only > 1:
                row["hap_classification"] = "germline"
            elif case_both > 1 and case_none > 1 and case_germ_only > 1 and case_var_only <= 1:
                row["hap_classification"] = "mosaic"
            elif case_both <= 1 and case_none > 1 and case_var_only > 1 and case_germ_only > 1:
                row["hap_classification"] = "mosaic"
            elif case_both <= 1 and case_var_only <= 1:
                row["hap_classification"] = "artifact"
            else:
                row["hap_classification"] = "artifact"

            results.append(row)

    for bam in bam_files:
        bam.close()
    return results

def write_phasing_tags(results, out_prefix):
    """Generate phasing_tags.tsv for bcftools annotate."""
    with open(f"{out_prefix}_phasing_tags.tsv", "w") as f:
        f.write("#CHROM\tPOS\tPHASING\n")
        for r in results:
            chrom = r["chrom"]
            pos = r["Var_pos"]
            phase = r["hap_classification"].upper()
            if phase == "MOSAIC":
                tag = "MOSAIC_PHASED"
            elif phase in ["GERMLINE", "ARTIFACT"]:
                tag = phase
            else:
                tag = "UNABLE_TO_PHASE"
            f.write(f"{chrom}\t{pos}\t{tag}\n")

def main():
    parser = argparse.ArgumentParser(description="Phasing step: read-level haplotype analysis, classification, and tag output.")
    parser.add_argument("-t", "--tsv", required=True, help="TSV from Step 4")
    parser.add_argument("-b", "--bams", nargs="+", required=True, help="PacBio BAMs")
    parser.add_argument("-s", "--sex", default="unknown", choices=["male", "female", "unknown"], help="Sex of the sample (for haploid X/Y logic)")
    parser.add_argument("-i", "--id", default="haplotype_summary", help="Output prefix for results")
    args = parser.parse_args()

    bam_list = args.bams
    results = count_haplotypes(args.tsv, bam_list, args.sex)

    # Write main classification table
    out_tsv = f"{args.id}.haplotyped.tsv"
    fieldnames = list(results[0].keys()) if results else []
    with open(out_tsv, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    # Write phasing tags for bcftools annotate
    write_phasing_tags(results, args.id)

    print(f"[Done] Classification TSV: {out_tsv}")
    print(f"[Done] Phasing tags TSV: {args.id}_phasing_tags.tsv")

if __name__ == "__main__":
    main()
