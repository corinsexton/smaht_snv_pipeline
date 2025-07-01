#!/usr/bin/python

import sys
import gzip
f_vcf = sys.argv[1]
f_out = sys.argv[2]
f = gzip.open(f_vcf, 'rt')
out = open(f_out, 'w')

cnt_all, cnt_rm, cnt_pass = 0, 0, 0

for line in f:
    if line[0] == '#':
        out.write(line)
        if line[:2] == '#C':
            header = line.strip().split('\t')
            print(header)
    else:
        cnt_all += 1
        s = line.strip().split('\t')
        l_variant = [s[0], s[1], s[3], s[4]]
        variant = ':'.join(l_variant)
        INFO = s[7].split(';')
        FILTER = s[6]
        popaf = None
        for info in INFO:
            if 'AF_grpmax_joint' in info:
                popaf = float(info.split('=')[1])
            elif 'AF_joint' in info:
                if popaf is not None:
                    popaf_joint = float(info.split('=')[1])
                    if popaf_joint > popaf:
                        popaf = popaf_joint
                else:
                    popaf = float(info.split('=')[1])

        if popaf is None or popaf < 0.001:
            cnt_pass += 1
            out.write(line)
        else:
            cnt_rm += 1

print("gnomAD 0.001 All, Passed, Removed: ", cnt_all, cnt_pass, cnt_rm)
