from collections import Counter
from scipy.stats import binomtest
import sys
import gzip

ID=sys.argv[1]
tool =sys.argv[2]
LRvalid_file = sys.argv[3]
f_input = sys.argv[4]
final_output = sys.argv[5]

dic_sampledepth = {"ST001-1A":255, 'ST001-1D':485, 'ST002-1D': 378,'ST002-1G': 104, 'ST003-1Q':196, 'ST004-1Q':156}
dic_tool = {"0.MF": "MF", "0.MT":"MT", "2.RUFUS":"RF", "4.STK":"STK", "1.MTto":"MTto"}
tool_short = dic_tool[tool]
path="/home/yoh855/jiny/SMaHT/Tissue6/truthset/2.LR_phase/"+tool+"/phasing/"
outpath = "/home/yoh855/jiny/SMaHT/Tissue6/truthset/2.LR_phase/"+tool+'/output/'
input_file = f_input
sample_depth = dic_sampledepth[ID]
variant_out =open(final_output, 'w')

# Open and read the input file
with open(input_file, "r") as f:
    # Initialize variables
    header = []
    data = []

    # Process each line in the file
    for line in f:
        # Strip whitespace and split by tab
        s = line.strip().split('\t')

        # Handle header line
        if s[0] == 'chrom':
            header = s
            header.append("hap_classification")
            header.append("cnt_germline")
            header.append("cnt_candidate")
            header.append("cnt_wildtype")
            header.append("Final_classification")
            header.append("Binom_p_SR")
            header.append("Binom_p_LR")
            header.append("Decision")
            header.append("Tag")
            print("Header:", header)  # Debug print for header
        else:
            dic = dict(zip(header[:-9], s))
            #print (dic)
            if dic['Germ_pos'] == "NA":
                dic['hap_classification'] = 'Not_applicable'
            else:
                total_read = float(dic['total_read'])
                case_both = float(dic['case_both'])
                case_var_only = float(dic['case_var_only'])
                case_germ_only = float(dic['case_germ_only'])
                case_none = float(dic['case_none']) #In case the pileup can't be found because of local realignment of the file
                cnt_germline = case_both + case_germ_only
                cnt_var = case_both + case_var_only
                dic['cnt_germline'] = cnt_germline
                dic['cnt_var'] = cnt_var
         
                dic['depth_sample_hap']=dic_sampledepth[ID]/2
                if dic['chrom'] != 'chrX' and dic['chrom'] != 'chrY':
                    if case_both > 1 and case_none > 1 and case_var_only <=1 and case_germ_only <=1: #FORGIVE readcount =1, as it could be error
                        dic['hap_classification'] = 'germline' #hap=2
                    elif case_both <= 1  and case_none <=1 and case_var_only > 1 and case_germ_only > 1 :
                        dic['hap_classification'] = 'germline' #hap=2
                    elif case_both >1  and case_none > 1 and case_germ_only > 1 and case_var_only <= 1:
                        dic['hap_classification'] = 'mosaic' #hap=3'
                    elif case_both <=1 and case_none > 1 and case_var_only > 1 and case_germ_only > 1:
                        dic['hap_classification'] = 'mosaic' #hap=3'
                    elif case_both <= 1 and case_var_only <= 1:
                        dic['hap_classification'] = 'artifact' #hap=novar'
                    else:
                        dic['hap_classification'] =  'artifact' # 'hap>3'
                    
                    dic['expected_hap_count'] = total_read/dic['depth_sample_hap'] if dic['depth_sample_hap'] == 0 else 'NA'
                    dic['expected_vaf_per_hap'] = 1/dic['expected_hap_count'] if dic['depth_sample_hap'] == 0 else 'NA'
                else: #on sex chromosome phased with nearby homSNP
                    if case_both > 1 and case_none<=1 and case_germ_only > 1 and case_var_only <= 1:
                        dic['hap_classification'] = 'mosaic'
                    else:
                        dic['hap_classification'] = 'artifact'

                    
                    
                    
            data.append(dic)

# Group by 'chrom' and 'Var_pos'
grouped_data = {}
for entry in data:
    key = (entry['chrom'], entry['Var_pos'])
    grouped_data.setdefault(key, []).append(entry)

final_output = []
classification_summary = Counter()
Tag_summary = Counter()
dic_vafinfo = {}

# Read variant file and extract VAF information
variant_header = []
variant_lines = []

with gzip.open(LRvalid_file, "rt") as variant_file:
    for line in variant_file:
        if line.startswith("##"):
            variant_header.append(line)
        elif line.startswith("#C"):
            variant_header.append(line)
            l_header = line.strip().split('\t')
            print(l_header)
        else:
            variant_lines.append(line)
            dic = dict(zip(l_header, line.strip().split('\t')))
            #SR_vaf, SR_dp, SR_alt, LR_vaf, LR_dp, LR_alt = dic['SR_vaf'], dic['SR_dp'], dic['SR_alt'], dic['LR_vaf'], dic['LR_dp'], dic['LR_alt']
            #l_vafinfo = [SR_vaf, SR_dp, SR_alt, LR_vaf, LR_dp, LR_alt]
            LR_vaf, LR_dp, LR_alt = dic['LR_vaf'], dic['LR_dp'], dic['LR_alt']
            l_vafinfo = [ LR_vaf, LR_dp, LR_alt]
            key = (dic['#CHROM'], dic['POS'])
            dic_vafinfo[key] = l_vafinfo

#print (data)
for entries in data:
    hap_classification = entries['hap_classification']
    classification_summary[hap_classification] += 1
    key = (entries['chrom'], entries['Var_pos'])
    vaf_info = dic_vafinfo.get(key, ["NA"] * 6)
    entries['Binom_relativeLR'] = None
    Tag = []
    try:
        #SR_dp = int(vaf_info[1])
        #SR_alt = int(vaf_info[2])
        #Binom_p_SR = binomtest(SR_alt, SR_dp, p=0.5).pvalue if SR_dp > 0 else "NA"
        #entries['Binom_p_SR'] = Binom_p_SR
        Binom_p_SR = "NA"
    except ValueError:
        Binom_p_SR = "NA"
        entries['Binom_p_SR'] = Binom_p_SR

    try:
        #print (vaf_info)
        LR_dp = int(vaf_info[1])
        LR_alt = int(vaf_info[2])
        Binom_p_LR = binomtest(LR_alt, LR_dp, p=0.5).pvalue if LR_dp > 0 else "NA"
        vaf_LR = float(vaf_info[0])
        entries['Binom_p_LR'] = Binom_p_LR
        
    except ValueError:
        Binom_p_LR = "NA"
        vaf_LR = "NA"
        entries['Binom_p_LR'] = Binom_p_LR

    if entries['hap_classification'] not in ['mosaic']:
        decision = 'Failed'
        Tag =['Phasing_fail']
        entries['Decision'] = decision

    else:
        decision = 'PASS'
        entries['Decision'] = decision
        Tag = []
        # (key, entries['total_read'], entries['depth_sample_hap'])
        if Binom_p_LR != "NA":
            if float(entries['Binom_p_LR']) > 0.01 or  float(LR_alt)/int(LR_dp) > 0.5:
                Tag.append("VAF_high")
        if entries['chrom'] != 'chrX' and entries['chrom'] != 'chrY':
            if float(entries['total_read']) < entries['depth_sample_hap']:
                Tag.append("Weak_align")
            if int(entries['total_read']) > 0:
                expected_vaf_for_onehap = entries['depth_sample_hap']/float(entries['total_read'])
                if int(entries['total_read'])>0 and expected_vaf_for_onehap > 1:
                    pass
                else:
                    binom_relative = binomtest(int(entries['cnt_var']), int(entries['total_read']), expected_vaf_for_onehap).pvalue 
                    entries['Binom_relativeLR'] = binom_relative
                    if binom_relative >= 0.01:
                        Tag.append("pGerm_segdup")
        else:
            pass
        if Tag == []:
            Tag = ['HighConf']

    entries['Tag'] = ';'.join(Tag)
    Tag = entries['Tag']
    Tag_summary[Tag] += 1

    final_output.append({
        'chrom': entries['chrom'],
        'Var_pos': entries['Var_pos'],
        'Var_ref': entries['Var_ref'],
        'Var_alt': entries['Var_alt'],
        'hap_classification' : entries['hap_classification'],
        'SR_vaf': '.', #vaf_info[0],
        'SR_dp': '.', #vaf_info[1],
        'SR_alt': '.', #vaf_info[2],
        'LR_vaf': vaf_info[0],
        'LR_dp': vaf_info[1],
        'LR_alt': vaf_info[2],
        'Binom_p_SR': 'NA', #entries['Binom_p_SR'],
        'Binom_p_LR':entries['Binom_p_LR'],
        'Binom_relativeLR' : entries['Binom_relativeLR'],
        'Decision': entries['Decision'],
        'Tag' : entries['Tag']
    })
dic_decision = { (entries['chrom'], entries['Var_pos']): entries['Decision'] for entries in final_output }
dic_Tag = { (entries['chrom'], entries['Var_pos']): entries['Tag'] for entries in final_output }
# Write filtered variants to output file
print (dic_Tag)
for header_line in variant_header:
    variant_out.write(header_line)
for line in variant_lines:
    dic = dict(zip(l_header, line.strip().split('\t')))
    key = (dic['#CHROM'], dic['POS'])
    if key in dic_decision and dic_Tag[key] == 'HighConf':
        print (dic_Tag[key])
        variant_out.write(line)

print(f"Filtered variant output written to {variant_out}")
'''
# Write the f#ull output to a new TSV file
output_file = input_file.replace(".tsv", "_classified.tsv")
with open(output_file, "w") as out:
    # Write header
    out.write("\t".join(header) + "\n")

    # Write data rows
    for entry in data:
        row = [str(entry[col]) for col in header]
        out.write("\t".join(row) + "\n")

print(f"Full output written to {output_file}")
'''
# Write the final output with selected columns
final_output_file = input_file.replace(".tsv", "_classified_final.tsv")
#print (final_output)
with open(final_output_file, "w") as out:
    # Write header
    final_header = ["chrom", "Var_pos", "Var_ref", "Var_alt",  "hap_classification","SR_vaf", "SR_dp", "SR_alt", "LR_vaf", "LR_dp", "LR_alt", "Binom_p_SR", "Binom_p_LR", 'Binom_relativeLR', "Decision", "Tag"]
    out.write("\t".join(final_header) + "\n")

    # Write final data rows
    for entry in final_output:
        row = [str(entry[col]) for col in final_header]
        out.write("\t".join(row) + "\n")

print(f"Final output written to {final_output_file}")

# Print classification summary
print("\nClassification Summary:")
for classification, count in classification_summary.items():
    print(f"{classification}: {count}")

print(f"Total entries: {len(final_output)}")

print("\nTag Summary:")
for Tag, count in Tag_summary.items():
    print(f"{Tag}: {count}")





