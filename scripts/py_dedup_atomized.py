#!/home/yoh855//miniconda3/bin/python

import pysam
import argparse

def filter_atomized_records(original_vcf, atomized_vcf, output_vcf):
    # Open the original and atomized VCF files (bgzip-compressed VCFs)
    original = pysam.VariantFile(original_vcf)  # Reading .vcf.gz
    atomized = pysam.VariantFile(atomized_vcf)  # Reading .vcf.gz
    
    # Create an output bgzipped VCF file
    output = pysam.VariantFile(output_vcf, 'w', header=atomized.header)
    
    # Store original records in a dictionary (key: CHROM:POS:REF:ALT)
    original_records = {}
    for record in original:
        key = f"{record.chrom}:{record.pos}:{record.ref}:{','.join(record.alts)}"
        if key not in original_records:
            original_records[key] = record  # Only keep one copy of the original record

    # Write the atomized records to the output, skipping any duplicate writes
    written_records = set()  # To track records we've already written

    for record in atomized:
        key = f"{record.chrom}:{record.pos}:{record.ref}:{','.join(record.alts)}"
        
        if key in original_records and key not in written_records:
            # Write the original record only once
            output.write(original_records[key])
            written_records.add(key)
            del original_records[key]  # Remove after writing to avoid duplicates
        elif key not in written_records:
            # If it's a new atomized record, write it
            output.write(record)
            written_records.add(key)
    
    # Close the VCF files
    original.close()
    atomized.close()
    output.close()

    # Compress the output VCF with bgzip and index it
    pysam.tabix_index(output_vcf, preset='vcf', force=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter out atomized duplicates and keep original records, using bgzip-compressed VCF files.")
    parser.add_argument("original_vcf", help="Original bgzipped VCF file (.vcf.gz)")
    parser.add_argument("atomized_vcf", help="Atomized bgzipped VCF file (.vcf.gz)")
    parser.add_argument("output_vcf", help="Filtered output bgzipped VCF file (.vcf.gz)")

    args = parser.parse_args()

    filter_atomized_records(args.original_vcf, args.atomized_vcf, args.output_vcf)

