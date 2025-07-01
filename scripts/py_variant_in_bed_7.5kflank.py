#!/usr/bin/python
import sys

# Path to the reference genome dictionary file
dict_file = "/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.dict"

# Build a dictionary mapping chromosome name -> chromosome length
chrom_lengths = {}
with open(dict_file, 'r') as df:
    for line in df:
        if line.startswith("@SQ"):
            # Split the line into fields (split on whitespace)
            fields = line.strip().split()
            chrom = None
            length = None
            for field in fields:
                if field.startswith("SN:"):
                    chrom = field.replace("SN:", "")
                elif field.startswith("LN:"):
                    try:
                        length = int(field.replace("LN:", ""))
                    except ValueError:
                        length = None
            if chrom is not None and length is not None:
                chrom_lengths[chrom] = length
#print (chrom_lengths['chrX'])
# Input and output file paths provided as command-line arguments
infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as f, open(outfile, 'w') as out:
    for line in f:
        if line.startswith('#'):
            continue  # skip header lines
        s = line.strip().split('\t')
        # Assume the first column is the chromosome and the second column is the position.
        chrom = s[0]
        pos = int(s[1])
        
        # Compute start: if pos > 7500, subtract 7500; otherwise 0.
        start = pos - 5000 if pos > 5000 else 0
        
        # Compute tentative end
        end = pos + 5000 if pos + 5000 < chrom_lengths[chrom] else chrom_lengths[chrom]

        
        # Write the output (chrom, start, end) separated by tabs.
        out.write(f"{chrom}\t{start}\t{end}\n")
