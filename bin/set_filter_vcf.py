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

    if 'PASS' not in new_header.filters:
        new_header.add_line('##FILTER=<ID=PASS,Description="At least one cross-evidence label present (CrossTech, CrossCaller, CrossTissue, or CrossCore)">')

    if 'LowEvidence' not in new_header.filters:
        new_header.add_line('##FILTER=<ID=LowEvidence,Description="No cross-evidence: lacks CrossTech, CrossCaller, CrossTissue, and CrossCore">')

    if 'EvidenceScore' not in new_header.info:
        new_header.add_line('##INFO=<ID=EvidenceScore,Number=1,Type=Integer,Description="Count of cross-evidence labels present (CrossTech, CrossCaller, CrossTissue, CrossCore)">')

    return new_header


###############################################################################
# Main
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Assign PASS / LowEvidence FILTERs and EvidenceScore to variants"
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

        crossTech   = "CrossTech"   in rec.info
        crossTissue = "CrossTissue" in rec.info
        crossCaller = "CrossCaller" in rec.info
        crossCore   = "CrossCore"   in rec.info

        evidence_score = sum([crossTech, crossTissue, crossCaller, crossCore])

        # clone record so filters can be added freely
        new = clone_record(rec, vcf_out.header)
        if 'TIER' in new.info:
            del new.info['TIER']

        new.info['EvidenceScore'] = evidence_score

        new.filter.clear()
        if evidence_score == 0:
            new.filter.add("LowEvidence")
        else:
            new.filter.add("PASS")

        vcf_out.write(new)

    vcf_out.close()
    print(f"[set_filter_vcf] Wrote output: {out_path}")


###############################################################################
# Entrypoint
###############################################################################
if __name__ == "__main__":
    main()

