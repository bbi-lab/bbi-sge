#!/usr/bin/env python

import os
import sys
import subprocess
import argparse


# Taking variables that are declared in main.nf by argparser
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='parameters for nomenclature') 
	parser.add_argument('--path', help='directory')
	parser.add_argument('--output', help='ouputtest')
	args = parser.parse_args()

working_dir = args.path

for i in os.listdir(working_dir):
    if i.endswith(".fastq.gz"): 
        index_of_first_dot = i.find('.')
        index_of_first_under = i.find('_')

        index_of_extension = i.find('.fastq.gz')
        sample_name = i[:index_of_first_dot]
        sample_lib = i[:index_of_first_under]
        before_extension_name = i[:index_of_extension]
        
    else:
        pass

