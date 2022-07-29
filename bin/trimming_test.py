#!/usr/bin/env python

import gzip
import sys
import os
import argparse
from itertools import islice

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parameters for sam to edits') 
    parser.add_argument('--path', help='amplicon list')
    parser.add_argument('--output', help='ouputtest')
    args = parser.parse_args()

working_dir = args.path
output_dir = args.output

#outfile_path = output_dir+'/'+name+'.txt'
#out_file = open(sys.argv[4],'w')




out_file = open(output_dir,'w')
        
with gzip.open(working_dir, 'r') as my_reads:
    reads = 0
    reads_w_n = 0
    while True:
        my_fastq = list(islice(my_reads, 4))
        if not my_fastq:
            break
        else:
            reads+=1
            seq = my_fastq[1]
            if 'N' in seq:
                reads_w_n += 1
            elif 'N' not in seq:
                out_file.write(''.join(my_fastq))
