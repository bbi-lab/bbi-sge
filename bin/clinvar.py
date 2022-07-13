#!/usr/bin/env python

import gzip
import pandas as pd
import argparse
import time

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='parsing chromosome and postions information')
        #parser.add_argument('--chr', type=int, help='chromosome #')
        args = parser.parse_args()

#nrows = 10001
#delim_whitespace=True

data = pd.read_csv('/mnt/c/Users/shinj/source/repos/bbi-sge/bin/variant_summary.txt.gz', seq='\s+', header=1, compression='gzip', nrows=10)

print(data)
# parsing input to #Chrome , start position #, and end position #
'''
chrome = data[data['#Chrom'] == args.chr ]
print (chrome)

pos = data[data['Pos'].between(args.start, args.end)]
print (pos)

filtering = data[(data['#Chrom'] == args.chr) & (data['Pos'].between(args.start, args.end))]
print (filtering)

filtering.to_csv('./test.csv')
'''
