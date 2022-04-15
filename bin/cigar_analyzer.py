#!/usr/bin/env python

import os
import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parameters for sam to edits') 
    parser.add_argument('--path', nargs='+', help='amplicon list')
    parser.add_argument('--cigar', help='experiment list') 
    parser.add_argument('--output', help='ouputtest')
    args = parser.parse_args()

def build_cigar_dict(sam_file):
    line_count = 0
    cigar_counts = {}
    for line in sam_file:
        if line_count < 2:
            line_count += 1
            continue
        else:
            line_count += 1
            #pull out cigar string
            cigar_string = line.split('\t')[5]
            #print cigar_string
            if cigar_string in cigar_counts:
                cigar_counts[cigar_string]+=1
            else:
                cigar_counts[cigar_string]=1
    return cigar_counts #returns the final dictionary after parsing the file

def sort_counts_dict(counts_dict):
    return sorted(counts_dict.keys(), key=lambda k: -1*int(counts_dict[k]))

def get_read_count(counts_dict):
    total_reads = 0
    for my_key in counts_dict:
        total_reads += counts_dict[my_key]
    return total_reads


def write_top_100_outfile(dict1,dict2,file_path):
    with open(file_path,'w') as f:
        f.write('cigar\ts1_reads\ts2_reads\ts1_rpt\ts2_rpt\ts2_s1_ratio\n')
        outlines = 0
        sorted_dict1_keys = sort_counts_dict(dict1)
        sorted_dict2_keys = sort_counts_dict(dict2)
        total_dict1_reads = get_read_count(dict1)
        total_dict2_reads = get_read_count(dict2)
        sorted_dict1_key_count = len(sorted_dict1_keys)
        max_lines = 100
        if sorted_dict1_key_count < 100:
            max_lines = sorted_dict1_key_count
        while outlines < max_lines:
            outlines += 1
            #calculate all those things above, make them strings, and write them out after tab-joining.
            my_key = sorted_dict1_keys[outlines-1]
            d5_reads = dict1[my_key]
            if my_key in dict2:
                d11_reads = dict2[my_key]
            else:
                d11_reads = 0
            d5_rpt = float(d5_reads)/total_dict1_reads*1000
            d11_rpt = float(d11_reads)/total_dict2_reads*1000
            if d5_rpt != 0:
                d11_d5_ratio = d11_rpt/d5_rpt
            else:
                d11_d5_ratio = 9999
            output_list = [str(my_key),str(d5_reads),str(d11_reads),str(d5_rpt),str(d11_rpt),str(d11_d5_ratio)]
            f.write('\t'.join(output_list)+'\n')


working_dir = args.path
#print working_dir
output_dir = args.output
#print output_dir
counts_dod = {}
#for i in working_dir:
for files in working_dir:
    print working_dir
    if files.endswith('.sam'):
        print files
        counts_dod[files] = build_cigar_dict(open(files,'r'))

comparison_list = args.cigar 
#print comparison_list
comps = comparison_list.split(',')
#print comps

for comp in comps:
    samples = comp.split('+')
    print samples 
    sample_1 = samples[0]
    print samples[0]
    sample_2 = samples[1]
    print sample_2
    write_top_100_outfile(counts_dod[sample_1+'.sam'],counts_dod[sample_2+'.sam'],output_dir+'/'+sample_1+'_'+sample_2+'_top100.txt')
