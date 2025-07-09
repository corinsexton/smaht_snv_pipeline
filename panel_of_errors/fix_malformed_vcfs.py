#!/usr/bin/env python

import sys

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Trim VCF-like columns based on colon-containing fields.")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input file (e.g. TSV or VCF)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output file where trimmed results will be saved."
    )
    return parser.parse_args()

args = parse_args()

full_header = ''
first_line = True


file = args.input
outfile = open(args.output,'w')

with open(file,'r') as infile:
    for line in infile:
        if line[0:2] == '##':
            full_header += line
            continue
        if line[0] == '#':
            outfile.write(full_header)
            header = line.strip().split()
        elif first_line:
            first_line = False
            ll = line.strip().split()
            last_colon_index = max(i for i, v in enumerate(ll) if ':' in v)

            new_header = header[:last_colon_index + 1]
            ll_trimmed = ll[:last_colon_index + 1]
            outfile.write("\t".join(new_header) + '\n')
            outfile.write("\t".join(ll_trimmed) + '\n')
        else:
            ll = line.strip().split()
            ll_trimmed = ll[:last_colon_index + 1]
            outfile.write("\t".join(ll_trimmed) + '\n')

outfile.close()
