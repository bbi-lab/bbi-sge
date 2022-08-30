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
data = pd.read_csv(args.path, sep=',', header=0, low_memory=False)

##bcl2 = data[]
r1 = data[data['fastq_1']]
print (r1)

r2 = data[data['fastq_2']]
print (r2)

