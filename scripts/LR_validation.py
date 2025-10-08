#!/usr/bin/env python3
import sys
import os
import gzip
import pysam
from collections import defaultdict

def parse_pileup_counts(pileup_vcf_path):
    """
    Parse the multisample pileup VCF. Returns:
      - counts[sample][(chrom,pos,ref,alt)] = [ref_count, alt_count, sec_count]
      - sample_ids: dict mapping sample-column -> tissue ID (e.g. ST001-1A)
    """
    vcf = pysam.VariantFile(pileup_vcf_path)
    samples = list(vcf.header.samples)
    counts = {s: {} for s in samples}
    sample_ids = {}
    for s in samples:
        # map full column name to just the tissue ID (first part before '_')
        base = os.path.basename(s)
        sid = base.split('_')[0]  # e.g. 'ST001-1A'
        sample_ids[s] = sid

    for rec in vcf:
        key = (rec.chrom, rec.pos, rec.ref, rec.alts[0] if rec.alts else None)
        for s in samples:
            ad = rec.samples[s].get('AD')
            ref_c = alt_c = sec_c = 0
            if ad and isinstance(ad, (tuple, list)):
                ref_c = int(ad[0])
                alt_c = int(ad[1]) if len(ad) > 1 else 0
                for other in ad[2:]:
                    cnt = int(other)
                    if cnt > sec_c:
                        sec_c = cnt
            counts[s][key] = [ref_c, alt_c, sec_c]
    return counts, sample_ids


def main():
    if len(sys.argv) != 4:
        print("Usage: python annotate_with_LR_multi.py processed.vcf.gz pileup.vcf.gz out_filtered.vcf")
        sys.exit(1)

    proc_vcf, pileup_vcf, out_vcf = sys.argv[1:]

    # derive fail output filename
    prefix = os.path.splitext(out_vcf)[0]
    fail_vcf = f"{prefix}OUT_LR.vcf"

    # extract case and donor from processed VCF filename
    case = os.path.basename(proc_vcf).split('_')[1].split('.')[0]  # e.g. ST001-1A
    main_donor = case.split('-')[0]

    # parse pileup counts
    counts, sample_ids = parse_pileup_counts(pileup_vcf)
    all_samples = list(counts.keys())



    # CS EDIT
    # identify target and suffix samples
    #target_col = None
    #for s in all_samples:
    #    if sample_ids[s] == case:
    #        target_col = s
    #        break
    #if target_col is None:
    #    sys.stderr.write(f"Error: no pileup sample matches case '{case}'\n")
    #    sys.exit(1)
    #suffix_cols = [s for s in all_samples if s != target_col]

    suffix_cols = [s for s in all_samples]

    # prepare output writers and write headers
    with gzip.open(proc_vcf, 'rt') as fin, \
         open(out_vcf, 'w') as pass_fout, \
         open(fail_vcf, 'w') as fail_fout:

        for line in fin:
            if line.startswith('##'):
                pass_fout.write(line)
                fail_fout.write(line)
            elif line.startswith('#'):
                #cols = ['LR_vaf','LR_dp','LR_alt','LR_sec_alt']
                cols = [f"{sample_ids[s]}_LRvaf" for s in suffix_cols]
                header = line.rstrip('\n') + '\t' + '\t'.join(cols) + '\n'
                pass_fout.write(header)
                fail_fout.write(header)
                break

    # iterate and annotate
    proc = pysam.VariantFile(proc_vcf)
    kept = failed = total = 0
    with open(out_vcf, 'a') as pass_fout, open(fail_vcf, 'a') as fail_fout:
        for rec in proc:
            total += 1
            key = (rec.chrom, rec.pos, rec.ref, rec.alts[0] if rec.alts else None)

            # CS EDIT
            ## target counts
            #r0, a0, s0 = counts[target_col].get(key, [0,0,0])
            ## check primary alt > secondary and >1
            #if not (a0 > s0 and a0 > 1):
            #    continue

            coverage = False
            overlap = False
            suffix_vafs = []
            for s in suffix_cols:
                r1, a1, _ = counts[s].get(key, [0,0,0])
                if (r1 + a1) > 0:
                    coverage = True
                    vaf1 = a1 / (r1 + a1)
                    suffix_vafs.append(f"{vaf1:.4f}")
                else:
                    suffix_vafs.append('.')
                donor = sample_ids[s].split('-')[0]
                if donor != main_donor and a1 > 1:
                    overlap = True

            ## require at least one other sample coverage
            #if not coverage:
            #    continue

            # CS EDIT
            #dp0 = r0 + a0
            #vaf0 = a0 / dp0 if dp0 else 0.0
            #sec_alt = s0

            line = str(rec).split('\n',1)[0]
            out_line = (
                f"{line}\t"
                + ''.join(f"\t{v}" for v in suffix_vafs)
                + "\n"
            )

            if overlap:
                fail_fout.write(out_line)
                failed += 1
            else:
                pass_fout.write(out_line)
                kept += 1

    sys.stdout.write(f"Total variants processed: {total}\n")
    sys.stdout.write(f"Kept (no overlap)       : {kept}\n")
    sys.stdout.write(f"Failed (overlap)        : {failed}\n")

if __name__ == '__main__':
    main()
