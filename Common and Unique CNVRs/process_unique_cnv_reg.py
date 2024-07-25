from pybedtools import BedTool
import pandas as pd
import re

def extract_unique_reg(x):
    global unique_cnv_bed_all
    print(x)
    target_cnv_reg = x['CNVID']
    target_cnv_reg_bed = re.sub(r'[:-]', '    ', target_cnv_reg)
    target_cnv_reg_bed = BedTool(target_cnv_reg_bed, from_string=True)
    common_reg_list = x['all_common_reg'].split(';')
    common_reg_list_bed = '\n'.join([re.sub(r'[:-]', '    ', i) for i in common_reg_list])
    common_reg_list_bed = BedTool(common_reg_list_bed, from_string=True)
    unique_cnv_bed = target_cnv_reg_bed.subtract(common_reg_list_bed)
    unique_cnv_bed_all = unique_cnv_bed_all.cat(unique_cnv_bed)
    unique_cnv_reg = [i.chrom + ':' + str(i.start) + '-' + str(i.end) for i in unique_cnv_bed]
    unique_cnv_reg = ';'.join(unique_cnv_reg)
    return unique_cnv_reg

def extract_unique_length(x):
    print(x)
    target_cnv_reg = x['CNVID']
    target_cnv_reg_bed = re.sub(r'[:-]', '    ', target_cnv_reg)
    target_cnv_reg_bed = BedTool(target_cnv_reg_bed, from_string=True)
    common_reg_list = x['all_common_reg'].split(';')
    common_reg_list_bed = '\n'.join([re.sub(r'[:-]', '    ', i) for i in common_reg_list])
    common_reg_list_bed = BedTool(common_reg_list_bed, from_string=True)
    unique_cnv_bed = target_cnv_reg_bed.subtract(common_reg_list_bed)
    
    unique_cnv_length = [(i.end-i.start+1) for i in unique_cnv_bed]
    unique_cnv_length = sum(unique_cnv_length)
    return unique_cnv_length

unique_cnv_bed_all = BedTool("chr1  1   10000", from_string=True)
common_reg = pd.read_excel("/Users/jadeliang/Desktop/Mona_project_CNV/cnv_bed_files/common_CNV_reg.xlsx")
common_reg_only = common_reg[["common_with_ASD", "common_with_ASD_ID_DD", "common_with_ID_DD", "common_with_SCZ"]]
common_reg['all_common_reg'] = common_reg_only.apply(lambda x: ';'.join(x.dropna().values.tolist()), axis=1)
all_common_reg = common_reg[['CNVID', 'all_common_reg']]
all_common_reg['unique_reg'] = all_common_reg.apply(lambda x: extract_unique_reg(x), axis=1)
all_common_reg['unique_sum_length'] = all_common_reg.apply(lambda x: extract_unique_length(x), axis=1)
print(all_common_reg)

all_common_reg.to_csv('/Users/jadeliang/Desktop/Mona_project_CNV/cnv_bed_files/all_common_reg.csv')
unique_cnv_bed_all.saveas("/Users/jadeliang/Desktop/Mona_project_CNV/cnv_bed_files/unique_cnv_bed.bed")