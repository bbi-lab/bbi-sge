#2_seqprep_BRCA1_pipeline.py
#to be run in the folder for the experiment with all the reads split by sample as fastq.gz files
#generates a shell script to run
#requires no underscores to be in sample names (later script will switch from dashes to underscores)
#outputs the stats from each seqprep call to a file called seqprep_stats.txt

import os
import sys
import subprocess
import argparse

work_dir = os.getcwd()

command_file_name = work_dir+"/run_seqprep.sh"
command_file = open(command_file_name, 'w')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parameters for sam to edits') 
    parser.add_argument('--path', nargs='+', help='amplicon list')
    parser.add_argument('--output', help='ouputtest')
    args = parser.parse_args()

working_dir = args.path

for i in working_dir():
    if i.endswith("R1_001.fastq.gz"): 
        index_of_first_under = i.find('_')
        index_of_R1 = i.find('R1')
        sample_name = i[:index_of_first_under]
        file_name = i[:index_of_R1]
        # here are adapters for F and R seq primers
        command_file.write(r'/net/gs/vol1/home/gf2/bin/SeqPrep/./SeqPrep -f '+
                           #file name... PALB2X11-lib_(S16)
                           file_name+r'R1_001.fastq.gz -r '+
                           file_name+r'R2_001.fastq.gz -1 Seqprep/R1/'+
                           #sample name...PALB2X11-lib
                           sample_name+r'.R1.fastq.gz -2 Seqprep/R2/'+
                           sample_name+r'.R2.fastq.gz -A GGTTTGGAGCGAGATTGATAAAGT -B CTGAGCTCTCTCACAGCCATTTAG -M 0.1 -s Seqprep/merged/'+
                           sample_name+'.merged.fastq.gz -m 0.001 -q 20 -o 20 &\n')
    else:
        pass
command_file.write("wait\n")
command_file.close()

