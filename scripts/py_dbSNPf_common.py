#!/usr/bin/python

import sys
import gzip
f_vcf = sys.argv[1]
f_out = sys.argv[2]
f = gzip.open(f_vcf, 'rt')
out = open(f_out, 'w')

cnt_all, cnt_rm, cnt_pass, cnt_none = 0, 0, 0, 0

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
                if ',' in info:
                    popaf = float(info.split('=')[1].split(',')[1])
                else:
                    popaf = float(info.split('=')[1])
            elif 'AF_joint' in info:
                if popaf is not None:
                    if ',' in info:
                        popaf_joint = float(info.split('=')[1].split(',')[1])
                        if popaf_joint > popaf:
                            popaf = popaf_joint
                    else:
                        popaf_joint = float(info.split('=')[1])
                        if popaf_joint > popaf:
                            popaf = popaf_joint
                else:
                    if ',' in info:
                        popaf = float(info.split('=')[1].split(',')[1])
                    else:
                        popaf = float(info.split('=')[1])
            elif 'TOPMED' in info:
                if ',' in info:
                    topmed_af = info.split('=')[1].split(',')[1]
                    if topmed_af != '.':
                        topmed = float(topmed_af)
                        #print (topmed_af)
                        if popaf != None:
                            if topmed > popaf:
                                popaf = topmed
                        else:
                            popaf = topmed
                else:
                    topmed_af = info.split('=')[1]
                    if topmed_af != '.':
                        topmed = float(topmed_af)
                        if popaf != None:
                            if topmed > popaf:
                                popaf = topmed
                        else:
                            popaf = topmed
                        
        
        if popaf != None:
            if popaf >= 0.001:
                cnt_pass += 1
                out.write(line)
            else:
                cnt_rm += 1
        else:
            cnt_none +=1 
            

print("dbSNP common >  0.001 All, Passed, Removed: ", cnt_all, cnt_pass, cnt_rm, cnt_none)
