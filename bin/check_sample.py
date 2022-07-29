#!/usr/bin/env python

import pandas as pd
import argparse
import time

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='parsing chromosome and postions information')
        parser.add_argument('--path', type=int, help='chromosome #')
        parser.add_argument('--out', type=int, help='start posiition #')
        args = parser.parse_args()

## 'Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'Sample_Project', 'Description', 'index', 'index2', 'GenomeFolder'
data = pd.read_csv(args.path, sep='\t', header=0, low_memory=False)

chrome = data[data['Sample_ID'] == args.chr ]
print ("testing chromosome")
print (chrome)

'''
# parsing input of Pos
pos = data[data['Pos'].between(args.start, args.end)]
print ("testing positions")
print (pos)

filtering = data[(data['#Chrom'] == args.chr) & (data['Pos'].between(args.start, args.end))]
print (filtering)

filtering.to_csv('./test.csv')
'''