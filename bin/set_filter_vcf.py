#!/usr/bin/env python3
import argparse
import pysam

###############################################################################
# Clone a pysam VariantRecord into a new record tied to output header
###############################################################################
def clone_record(rec, out_header):
    """
    Duplicate a record into the output header.
    Copies all INFO and FORMAT/sample fields; FILTER is cleared (overwritten by caller).
    Handles both single-sample and multi-sample VCFs.
    """
    new_rec = out_header.new_record(
        contig=rec.contig,
        start=rec.start,
        stop=rec.stop,
        id=rec.id,
        qual=rec.qual,
        alleles=rec.alleles,
        filter=None,
    )

    # Copy all INFO fields present in the output header
    for key, val in rec.info.items():
        try:
            new_rec.info[key] = val
        except Exception:
            # pysam reads Number=1 String fields containing commas as a tuple;
            # rejoin with commas to reconstruct the original string and retry.
            if isinstance(val, tuple):
                try:
                    new_rec.info[key] = ','.join(str(v) for v in val)
                except Exception:
                    pass

    # Copy FORMAT/sample fields for every sample, field by field so one
    # problematic field doesn't silently drop all FORMAT data for a sample
    for sample in rec.samples:
        for key, val in rec.samples[sample].items():
            try:
                new_rec.samples[sample][key] = val
            except Exception:
                pass

    return new_rec


def fix_header(header):
    """
    Return a copy of the input header with FILTER definitions added.
    All FORMAT, INFO, contig, and sample records from the input are preserved.
    """
    new_header = header.copy()

    filter_defs = [
        ('HighConf',      'High confidence variant (CrossTech is set, or any core has CrossCaller and CrossTissue is set)'),
        ('LowConf',       'Low confidence variant (any core has CrossCaller, or CrossTissue is set)'),
        ('LikelyArtifact','Lowest confidence: no CrossTech, CrossCaller, or CrossTissue evidence'),
    ]

    for flt_id, flt_desc in filter_defs:
        if flt_id not in new_header.filters:
            new_header.add_line(f'##FILTER=<ID={flt_id},Description="{flt_desc}">')

    return new_header


###############################################################################
# Main
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Assign HighConf / LowConf / . FILTERs to variants"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input VCF (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output VCF (.vcf or .vcf.gz). Default: <input>.confidence.vcf.gz"
    )
    args = parser.parse_args()

    in_path = args.input

    # Determine default output path
    if args.output:
        out_path = args.output
    else:
        if in_path.endswith(".vcf.gz"):
            out_path = in_path.replace(".vcf.gz", ".confidence.vcf.gz")
        elif in_path.endswith(".vcf"):
            out_path = in_path.replace(".vcf", ".confidence.vcf")
        else:
            out_path = in_path + ".confidence.vcf.gz"

    ###########################################################################
    # Open input VCF and modify header to add filter definitions
    ###########################################################################
    vcf_in = pysam.VariantFile(in_path)
    header = vcf_in.header.copy()

    new_header = fix_header(header)

    ###########################################################################
    # Create output VCF with updated header
    ###########################################################################
    vcf_out = pysam.VariantFile(out_path, "w", header=new_header)

    ###########################################################################
    # Process variants
    ###########################################################################
    for rec in vcf_in:

        # Tissue-level INFO flags
        crossTech   = "CrossTech"   in rec.info
        crossTissue = "CrossTissue" in rec.info

        # CrossCaller is an INFO flag (2+ unique callers across all cores)
        any_cross_caller = "CrossCaller" in rec.info

        # clone record so filters can be added freely
        new = clone_record(rec, vcf_out.header)
        if 'TIER' in new.info:
            del new.info['TIER']

        # Decision tree (plan section 3g):
        # HighConf     — CrossTech is set OR (any core has CrossCaller AND CrossTissue is set)
        # LowConf      — any core has CrossCaller OR CrossTissue is set
        # LikelyArtifact — none of the above

        new.filter.clear()

        if crossTech or (any_cross_caller and crossTissue):
            new.filter.add("HighConf")
        elif any_cross_caller or crossTissue:
            new.filter.add("LowConf")
        else:
            new.filter.add("LikelyArtifact")

        vcf_out.write(new)

    vcf_out.close()
    print(f"[set_filter_vcf] Wrote output: {out_path}")


###############################################################################
# Entrypoint
###############################################################################
if __name__ == "__main__":
    main()

