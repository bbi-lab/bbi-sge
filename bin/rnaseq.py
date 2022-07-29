#!/usr/bin/env python


import argparse
import sys
from itertools import islice

#parser arguments that are called in main.nf
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parameters for sam to edits') 
    parser.add_argument('--input', help='input files from trimming process')
    parser.add_argument('--output', help='ouput files for needleall')
    parser.add_argument('--amp', help='amp')
    args = parser.parse_args()




out_file = open(args.output,'w')
#summary_out = open('summary_cDNA_to_gDNA.txt','a')
cDNA_info_file = open('/net/bbi/vol1/home/jongs2/SGE/PALB2/PALB2_cDNA_data.txt','r') #tab-delimited file with header

#read in the editing data for all amplicons
line_count = 0
cDNA_info_dict = {} #  amp --> []
for line in cDNA_info_file: #file formated to be tdl:  exa.     amp    5_cDNA   3_cDNA  length  5_gDNA  3_gDNA  
    if line_count == 0:
        line_count +=1
    else:
        cDNA_info = line.strip().split()
        cDNA_info_dict[cDNA_info[0]] = cDNA_info[1:] #amp points to all other things


with open(sys.argv[2], 'r') as my_reads:
    amp = sys.argv[6]
    cDNA_5 = cDNA_info_dict[amp][0]
    cDNA_3 = cDNA_info_dict[amp][1]
    cDNA_len = int(cDNA_info_dict[amp][2])
    mod_5 = cDNA_info_dict[amp][3]
    mod_3 = cDNA_info_dict[amp][4]
    new_qs_5 = len(mod_5)*'H'
    new_qs_3 = len(mod_3)*'H'
    reads = 0
    full_length_reads = 0
    correct_5_reads = 0
    correct_3_reads = 0
    correct_5_3_reads = 0
    correct_5_3_full_length_reads = 0 
    while True:
        my_fastq = list(islice(my_reads, 4))
        if not my_fastq:
            break
        else:
            reads+=1
            seq = my_fastq[1].strip()
            qs = my_fastq[3].strip()
            if seq.startswith(cDNA_5):
                correct_5_reads +=1
                correct_5 = True
            else:
                correct_5 = False
            if seq.endswith(cDNA_3):
                correct_3_reads += 1
                correct_3 = True
            else:
                correct_3 = False
            if len(seq) == cDNA_len:
                full_length_reads += 1
                full_length = True
            else:
                full_length = False
                #print len(seq), 'length seq', cDNA_len
            if correct_5 and correct_3:
                correct_5_3_reads +=1
                correct_5_3 = True
            else:
                correct_5_3 = False
            if correct_5_3 and full_length:
                correct_5_3_full_length_reads +=1
                new_seq = mod_5+seq[len(cDNA_5):-len(cDNA_3)]+mod_3
                new_qs = new_qs_5+qs[len(cDNA_5):-len(cDNA_3)]+new_qs_3
                out_file.write(''.join([my_fastq[0],new_seq,'\n',my_fastq[2],new_qs,'\n']))
            else:
                pass
print sys.argv[2], 'finished processing.'
summary_out.write(sys.argv[1]+'\n'+str(reads)+' total reads processed.\n'+str(correct_5_3_full_length_reads)+" reads with correct 5' and 3' ends.\n"+str(correct_5_3_full_length_reads/float(reads))+' retention fraction.\n')
print reads, 'total reads processed'
print full_length_reads, 'full_length_reads'
print correct_5_reads, 'correct_5_reads'
print correct_3_reads, 'correct_3_reads'
print correct_5_3_reads, 'correct_5_3_reads'
print correct_5_3_full_length_reads, 'correct_5_3_full_length_reads (total output reads)'
print correct_5_3_full_length_reads/float(reads), 'retention fraction'











