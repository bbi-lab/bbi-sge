#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
import subprocess

# Taking variables that are declared in main.nf by argparser
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='parameters for nomenclature') 
	parser.add_argument('--path', help='sample sheet')
	parser.add_argument('--output', help='ouput')
	args = parser.parse_args()

out_file = open(sys.argv[4],'w')


with open(csv_file, 'rb') as csvfile:

    # get number of columns
    for line in csvfile.readlines():
        array = line.split(',')
        first_item = array[0]

    num_columns = len(array)
    csvfile.seek(0)

    reader = csv.reader(csvfile, delimiter=' ')
        included_cols = [1, 2, 6, 7]

    for row in reader:
            content = list(row[i] for i in included_cols)
            print content

'''
with open(sys.argv[2], 'r') as my_reads:
    fields = ['star_name', 'ra']

    df = pd.read_csv(args., skipinitialspace=True, usecols=fields)
    df = pd.read_csv("student.csv", usecols = ['Hours','Scores','Pass'])
print(df)
# See the keys
print df.keys()
# See content in 'star_name'
print df.star_name
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
    print sys.argv[1], 'finished processing.'
    print reads_w_n, 'reads with N base'
    print reads, 'total reads'
'''