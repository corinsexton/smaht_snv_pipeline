#!/usr/bin/env python3
import argparse
import pysam

###############################################################################
# Clone a pysam VariantRecord into a new record tied to output header
###############################################################################
def clone_record(rec, out_header):
    """
    Create a new record in the output VCF that duplicates the input record.
    This ensures FILTER + INFO fields added to the output header are valid.
    """

    key_order = [
        "CrossTech",
        "CrossCaller",
        "CrossTissue",
        "CALLERS",
        "ORIGINAL_FILTER",
        "CLUSTER",
        "CLUSTER_N",
        "PB_READ_CUTOFF",
        "SR_READ_CUTOFF",
        "ALT_SUPPORT",
        "SR_VAF",
        "TISSUE_PB_VAF",
        "TISSUE_ONT_VAF",
        "POOLED_PB_VAF",
        "POOLED_ONT_VAF",
        "TISSUE_SR_VAFS",
        "SR_ADF",
        "SR_ADR",
        "PB_ADF",
        "PB_ADR",
        "ONT_ADF",
        "ONT_ADR",
        "SB_SRC",
        "SB_PVAL",
        "FISHER",
        "GERMLINE_PVAL",
        "GERMLINE_PVAL_SR",
        "GERMLINE_PVAL_PB",
        "GERMLINE_PVAL_ONT",
        "GERMLINE_BINOM",
        "PB_PHASING",
        "REGION"
    ]

    info = dict(rec.info)
    ordered_info = {
        k: info[k]
        for k in key_order
        if k in info
    }


    new_rec = out_header.new_record(
        contig=rec.contig,
        start=rec.start,
        stop=rec.stop,
        id=rec.id,
        qual=rec.qual,
        alleles=rec.alleles,
        filter=None               # we overwrite FILTER
    )

    # Set INFO fields individually so pysam uses the correct Number/Type per field
    # (passing info= to new_record fails for Number=2 fields like SR_ADF)
    for k, v in ordered_info.items():
        new_rec.info[k] = v

    # copy FORMAT/sample fields
    for sample in rec.samples:
        new_rec.samples[sample].update(rec.samples[sample].items())

    return new_rec

def fix_header(header):
    """
    Create final output header
    """

    final_headers= [
            '##INFO=<ID=CrossTech,Number=0,Type=Flag,Description="Alt supported in both tissue short read and pooled PacBio data at or above their combined read thresholds">',
            '##INFO=<ID=CrossCaller,Number=0,Type=Flag,Description="Alt found in more than one variant caller">',
            '##INFO=<ID=CrossTissue,Number=0,Type=Flag,Description="Alt has VAF > 0 in another short read tissue">',
            '##INFO=<ID=CALLERS,Number=.,Type=String,Description="List of variant callers that reported this variant">',

            '##INFO=<ID=ORIGINAL_FILTER,Number=.,Type=String,Description="Original filter values">',
            '##INFO=<ID=CLUSTER,Number=1,Type=String,Description="Proximity clustering within window bp (PASS=not clustered, FAIL=clustered)">',
            '##INFO=<ID=CLUSTER_N,Number=1,Type=Integer,Description="If CLUSTER=FAIL, number of variants in the proximity cluster">',

            '##INFO=<ID=PB_READ_CUTOFF,Number=1,Type=Float,Description="Number of PacBio reads with ALT support required to pass">',
            '##INFO=<ID=SR_READ_CUTOFF,Number=1,Type=Float,Description="Number of Illumina reads with ALT support required to pass">',
            '##INFO=<ID=ALT_SUPPORT,Number=1,Type=String,Description="PASS/FAIL based on alt read support">',


            '##INFO=<ID=SR_VAF,Number=1,Type=Float,Description="VAF for short read in current tissue">',
            '##INFO=<ID=TISSUE_PB_VAF,Number=1,Type=Float,Description="VAF for PacBio in current tissue (if available)">',
            '##INFO=<ID=TISSUE_ONT_VAF,Number=1,Type=Float,Description="VAF for ONT in current tissue (if available)">',
            '##INFO=<ID=POOLED_PB_VAF,Number=1,Type=Float,Description="VAF for PacBio in current donor pooled tissues">',
            '##INFO=<ID=POOLED_ONT_VAF,Number=1,Type=Float,Description="VAF for ONT in current donor pooled tissues">',
            '##INFO=<ID=TISSUE_SR_VAFS,Number=.,Type=String,Description="VAFs for all tissues with short read nonzero VAF, based on pileup (BQ≥30)">',

            '##INFO=<ID=SR_ADF,Number=2,Type=Integer,Description="Tissue short read forward depths (REF,ALT)">',
            '##INFO=<ID=SR_ADR,Number=2,Type=Integer,Description="Tissue short read reverse depths (REF,ALT)">',
            '##INFO=<ID=PB_ADF,Number=2,Type=Integer,Description="Donor pooled Long-read forward depths (REF,ALT)">',
            '##INFO=<ID=PB_ADR,Number=2,Type=Integer,Description="Donor pooled Long-read reverse depths (REF,ALT)">',
            '##INFO=<ID=ONT_ADF,Number=2,Type=Integer,Description="Donor pooled ONT forward depths (REF,ALT)">',
            '##INFO=<ID=ONT_ADR,Number=2,Type=Integer,Description="Donor pooled ONT reverse depths (REF,ALT)">',

            '##INFO=<ID=SB_SRC,Number=1,Type=String,Description="Counts source used for Fisher strand test: PB, ONT, or SR">',
            '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher p-value for strand balance on chosen sample">',
            '##INFO=<ID=FISHER,Number=1,Type=String,Description="PASS/FAIL Fisher test">',
            '##INFO=<ID=GERMLINE_PVAL,Number=1,Type=Float,Description="Minimum binomial p-value for germline deviation across all platforms tested">',
            '##INFO=<ID=GERMLINE_PVAL_SR,Number=1,Type=Float,Description="Binomial p-value for germline deviation in tissue short read data">',
            '##INFO=<ID=GERMLINE_PVAL_PB,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled PacBio data">',
            '##INFO=<ID=GERMLINE_PVAL_ONT,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled ONT data">',
            '##INFO=<ID=GERMLINE_BINOM,Number=1,Type=String,Description="PASS/FAIL germline binomial test">',

            '##INFO=<ID=PB_PHASING,Number=1,Type=String,Description="Phasing classification from pooled PacBio nearest germline SNV haplotyping">',
            '##INFO=<ID=REGION,Number=1,Type=String,Description="SMaHT region classification: extreme, difficult, or easy">',

            '##FILTER=<ID=HighConf,Description="High confidence variant (CrossTech or CrossCaller+CrossTissue)">',
            '##FILTER=<ID=LowConf,Description="Low confidence variant (CrossCaller or CrossTissue only)">',
            '##FILTER=<ID=LikelyArtifact,Description="Variants passing all filters but with no CrossTech, CrossCaller, or CrossTissue evidence, lowest confidence variants">',
            '##FILTER=<ID=FAIL,Description="Variants failing one or more filters">'
    ]

    header_list = str(header).split('\n')

    new_header = pysam.VariantHeader()
    for header_line in header_list:
        if "fileformat" in header_line or "contig" in header_line:
            new_header.add_line(header_line)
        if "SAMPLE" in header_line:
            sample_line = header_line

    for header_line in final_headers:
        new_header.add_line(header_line)

    new_header.add_line(sample_line)

    return new_header


###############################################################################
# Main
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Assign HighConf / LowConf / LikelyArtifact / FAIL FILTERs to variants"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input VCF (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output VCF (.vcf or .vcf.gz). Default: <input>.confidence.vcf.gz"
    )
    parser.add_argument("--easy_regions", default=None, help="BED.gz for easy regions (tabix-indexed)")
    parser.add_argument("--diff_regions", default=None, help="BED.gz for difficult regions (tabix-indexed)")
    parser.add_argument("--ext_regions",  default=None, help="BED.gz for extreme regions (tabix-indexed)")
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

    # Open region tabix files if provided
    ext_tbx  = pysam.TabixFile(args.ext_regions)  if args.ext_regions  else None
    diff_tbx = pysam.TabixFile(args.diff_regions) if args.diff_regions else None
    easy_tbx = pysam.TabixFile(args.easy_regions) if args.easy_regions else None

    ###########################################################################
    # Create output VCF with updated header
    ###########################################################################
    vcf_out = pysam.VariantFile(out_path, "w", header=new_header)

    ###########################################################################
    # Process variants
    ###########################################################################
    for rec in vcf_in:

        # detection flags
        crossTech   = "CrossTech"   in rec.info
        crossCaller = "CrossCaller" in rec.info
        crossTissue = "CrossTissue" in rec.info

        # assign REGION (extreme > difficult > easy)
        region = None
        if ext_tbx:
            if any(True for _ in ext_tbx.fetch(rec.chrom, rec.pos, rec.pos + 1)):
                region = "extreme"
        if region is None and diff_tbx:
            if any(True for _ in diff_tbx.fetch(rec.chrom, rec.pos, rec.pos + 1)):
                region = "difficult"
        if region is None and easy_tbx:
            if any(True for _ in easy_tbx.fetch(rec.chrom, rec.pos, rec.pos + 1)):
                region = "easy"

        # clone record so filters can be added freely
        new = clone_record(rec, vcf_out.header)
        if region is not None:
            new.info["REGION"] = region

        # Decision tree:
        # If CrossTech OR (CrossCaller AND CrossTissue) → HighConf
        # Else if CrossCaller OR CrossTissue → LowConf
        # Else → LikelyArtifact

        new.filter.clear()

        info_string = str(rec).split('\t')[7]
        pb_phasing = rec.info.get("PB_PHASING", "")
        if 'FAIL' in info_string or pb_phasing in ("GERMLINE", "ARTIFACT","GERMLINE_SEGDUP"):
            new.filter.add("FAIL")
        else:
            if crossTech or (crossCaller and crossTissue):
                new.filter.add("HighConf")
            elif crossCaller or crossTissue:
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

