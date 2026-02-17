#!/usr/bin/env python3
import argparse
import pysam


def _ensure_info_fields(header: pysam.VariantHeader) -> None:
    """
    Add INFO fields if they don't already exist:
      - CLUSTER: PASS/FAIL
      - CLUSTER_N: number of variants in the proximity cluster (only set when FAIL)
    """
    if "CLUSTER" not in header.info:
        header.info.add(
            "CLUSTER",
            number=1,
            type="String",
            description="Proximity clustering within window bp (PASS=not clustered, FAIL=clustered)",
        )
    if "CLUSTER_N" not in header.info:
        header.info.add(
            "CLUSTER_N",
            number=1,
            type="Integer",
            description="If CLUSTER=FAIL, number of variants in the proximity cluster",
        )


def filter_clustered(input_vcf: str, output_vcf: str, window: int = 50) -> None:
    vcf_in = pysam.VariantFile(input_vcf)

    # Copy header and add our INFO fields
    header = vcf_in.header.copy()
    _ensure_info_fields(header)

    mode = "wz" if output_vcf.endswith(".gz") else "w"
    vcf_out = pysam.VariantFile(output_vcf, mode, header=header)

    group = []   # records in the current proximity group
    prev = None

    def flush_group(g):
        if not g:
            return
        clustered = (len(g) > 1)
        n = len(g)

        for rec_in in g:
            rec = vcf_out.new_record(
                contig=rec_in.contig,
                start=rec_in.start,
                stop=rec_in.stop,
                id=rec_in.id,
                alleles=rec_in.alleles,
                qual=rec_in.qual,
                filter=list(rec_in.filter.keys()) if rec_in.filter is not None else None,
                info=dict(rec_in.info),
            )

            rec.info["CLUSTER"] = "FAIL" if clustered else "PASS"
            if clustered:
                rec.info["CLUSTER_N"] = n
            else:
                # Only include CLUSTER_N on FAIL, per request
                if "CLUSTER_N" in rec.info:
                    del rec.info["CLUSTER_N"]
            vcf_out.write(rec)

    for rec in vcf_in.fetch():
        if prev is None or rec.chrom != prev.chrom or (rec.pos - prev.pos) > window:
            flush_group(group)
            group = [rec]
        else:
            group.append(rec)
        prev = rec

    flush_group(group)

    vcf_in.close()
    vcf_out.close()

    if output_vcf.endswith(".gz"):
        pysam.tabix_index(output_vcf, preset="vcf", force=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate proximity clusters within N bp via INFO fields (CLUSTER, CLUSTER_N)."
    )
    parser.add_argument("input_vcf", help="Input VCF file (bgzipped or not)")
    parser.add_argument("output_vcf", help="Output VCF file (bgzipped if ends with .gz)")
    parser.add_argument("--window", type=int, default=50, help="Window size in bp [default=50]")
    args = parser.parse_args()

    filter_clustered(args.input_vcf, args.output_vcf, args.window)

