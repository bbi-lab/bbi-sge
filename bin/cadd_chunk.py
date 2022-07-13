#!/usr/bin/env python

import gzip
import pandas as pd
import argparse
import itertools as it
import time

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='parsing chromosome and postions information') 
	parser.add_argument('--chr', type=int, help='chromosome #')
	parser.add_argument('--start', type=int, help='start posiition #')
	parser.add_argument('--end', type=int, help='end position #')
	args = parser.parse_args()

chunksizes = 10 ** 5
data = pd.read_csv('/net/bbi/vol1/cadd/nobackup/whole_genome_SNVs_inclAnno.tsv.gz', sep='\t', header=1, compression='gzip', chunksize=chunksizes, low_memory=False)
for chunk in data:
    print(chunk)
    
    #testing getting corresponding chromosome
    #chrome = chunk[chunk['#Chrom'] == args.chr ]
    #print ("testing for chromosome")
    #print (chrome)
    
    #testing getting corresponding position from start and end
    #pos = chunk[chunk['Pos'].between(args.start, args.end)]
    #print ("testing for positions")
    #print (pos)

    #filtering = chunk[(chunk['#Chrom'] == args.chr) & (chunk['Pos'].between(args.start, args.end))]
    #filtering.to_csv('./test.csv')
# parsing input to #Chrome , start position #, and end position #
'''
chrome = data[data['#Chrom'] == args.chr ]
print (chrome)

pos = data[data['Pos'].between(args.start, args.end)]
print (pos)

filtering = data[(data['#Chrom'] == args.chr) & (data['Pos'].between(args.start, args.end))]
print (filtering)
'''
