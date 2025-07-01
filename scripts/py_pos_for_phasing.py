#!/usr/bin/python

import sys

variant = sys.argv[1]
germline = sys.argv[2]

outname_raw = variant[:-4] + '.pos_for_phasing_raw.txt'
outname = variant[:-4] + '.pos_for_phasing.txt'

f_variant = open(variant, 'r')
f_germline = open(germline, 'r')
out_raw = open(outname_raw, 'w')
out = open(outname, 'w')

# Write headers
header = "chrom\tVar_pos\tVar_ref\tVar_alt\tGerm_pos\tGerm_ref\tGerm_alt"
out_raw.write(header + '\n')

dic_germline = {}
for line in f_germline:
    if line[0] != '#':
        s = line.strip().split('\t')
        chrom, pos, ref, alt = s[0], s[1], s[3], s[4]
        l_var = [chrom, pos, ref, alt]
        variant_id = ':'.join(l_var)
        flank = ','.join(s[-3:])
        if flank not in dic_germline:
            dic_germline[flank] = [variant_id]
        else:
            dic_germline[flank].append(variant_id)

for line in f_variant:
    if line[0] != '#':
        s = line.strip().split('\t')
        chrom, pos, ref, alt = s[0], s[1], s[3], s[4]
        l_var = [chrom, pos, ref, alt]
        variant_id = ':'.join(l_var)
        flank = ','.join(s[-3:])
        if flank in dic_germline:
            l_germline = dic_germline[flank]
            for germline_variant in l_germline:
                if germline_variant != variant_id:
                    # germline_variant is in the format "chrom:pos:ref:alt"
                    germ_info = germline_variant.split(':')[1:]  # use only pos, ref, alt
                    outline = [chrom, pos, ref, alt] + germ_info
                    out_raw.write('\t'.join(outline) + '\n')
        else:
            outline = [chrom, pos, ref, alt, 'NA', 'NA', 'NA']
            out_raw.write('\t'.join(outline) + '\n')

f_variant.close()
f_germline.close()
out_raw.close()

# Now, without using pandas, we group the raw output by (chrom, Var_pos)
# and for each group, select the row with the smallest absolute difference
# between Germ_pos and Var_pos.

# We'll use a dictionary: key = (chrom, Var_pos), value = (line, diff)
# If Germ_pos is "NA", we treat diff as None.
group_dict = {}

with open(outname_raw, 'r') as fin:
    header_line = fin.readline().strip()  # read header
    for line in fin:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        # fields: chrom, Var_pos, Var_ref, Var_alt, Germ_pos, Germ_ref, Germ_alt
        chrom, var_pos, germ_pos = fields[0], fields[1], fields[4]
        try:
            var_pos_int = int(var_pos)
        except:
            continue
        if germ_pos == "NA":
            diff = None
        else:
            try:
                germ_pos_int = int(germ_pos)
                diff = abs(germ_pos_int - var_pos_int)
            except:
                diff = None
        key = (chrom, var_pos)
        # If key is not yet stored, or if diff is valid and smaller than stored one, update.
        if key not in group_dict:
            group_dict[key] = (line, diff)
        else:
            stored_line, stored_diff = group_dict[key]
            # If current diff is not None and either stored_diff is None or current is smaller, update.
            if diff is not None:
                if stored_diff is None or diff < stored_diff:
                    group_dict[key] = (line, diff)
            # Otherwise, if diff is None, keep the stored value.
                    
# Write the filtered output file.
# Optionally, sort the keys for a predictable order.
sorted_keys = sorted(group_dict.keys(), key=lambda x: (x[0], int(x[1])))


out.write(header + "\n")
for key in sorted_keys:
    line, diff = group_dict[key]
    out.write(line + "\n")
