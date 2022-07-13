#!/usr/bin/env python

import pandas as pd
import argparse
import time

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='parsing chromosome and postions information') 
	parser.add_argument('--chr', type=int, help='chromosome #')
	parser.add_argument('--start', type=int, help='start posiition #')
	parser.add_argument('--end', type=int, help='end position #')
	args = parser.parse_args()

#nrows=10001

data = pd.read_csv('/net/bbi/vol1/cadd/nobackup/whole_genome_SNVs_inclAnno.tsv.gz', sep='\t', header=1, compression='gzip', nrows=10001, low_memory=False)
#skip rows?
#import csv - rather than pd
# trt ...

# parsing input to #Chrome , start position #, and end position 
# parsing input of #Chrom 
chrome = data[data['#Chrom'] == args.chr ]
print ("testing chromosome")
print (chrome)

# parsing input of Pos 
pos = data[data['Pos'].between(args.start, args.end)]
print ("testing positions")
print (pos)

filtering = data[(data['#Chrom'] == args.chr) & (data['Pos'].between(args.start, args.end))]
print (filtering)

filtering.to_csv('./test.csv')
