#!/usr/bin/env python3

"""
renaming with publishDir (nextflow)


working_dir = os.path.dirname(os.path.realpath(__file__))

def working_dir(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
"""

import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='removingi S after bcl2fastq') 
    parser.add_argument('--path', type=directory path)


for i in os.listdir(args.path):
    if i.endswith("fastq.gz"): #will only work on sam files
        index_of_first_S = i.find('S') ### change counts _ and S so that you can have S in name on the experiment
        index_of_first_under_post_S = i.find('_', index_of_first_S)
        sample_name = i[:index_of_first_S]+i[index_of_first_under_post_S+1:]
    else:
        pass

