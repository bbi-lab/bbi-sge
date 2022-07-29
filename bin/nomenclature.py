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

sample_mapping_dict = {}
with open(args.path, "r") as csv:

    ## Check header
    MIN_COLS = 2
    HEADER = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'Sample_Project', 'Description', 'index', 'index2', 'GenomeFolder']
    header = [x.strip('"') for x in csv.readline().strip().split(",")]
    if header[: len(HEADER)] != HEADER:
        print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
        sys.exit(1)

    ## Check sample entries
    for line in csv:
        lspl = [x.strip().strip('"') for x in line.strip().split(",")]

        # Check valid number of columns per row
        if len(lspl) < len(HEADER):
            print_error(
                "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                "Line",
                line,
            )
        num_cols = len([x for x in lspl if x])
        if num_cols < MIN_COLS:
            print_error(
                "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                "Line",
                line,
            )

        ## Check sample name entries
        sample, fastq_1, fastq_2 = lspl[: len(HEADER)]
        sample = sample.replace(" ", "_")
        if not sample:
            print_error("Sample entry has not been specified!", "Line", line)

        ## Check FastQ file extension
        for fastq in [fastq_1, fastq_2]:
            if fastq:
                if fastq.find(" ") != -1:
                    print_error("FastQ file contains spaces!", "Line", line)
                if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                    print_error(
                        "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                        "Line",
                        line,
                    )

        ## Auto-detect paired-end/single-end
        sample_info = []  ## [single_end, fastq_1, fastq_2]
        if sample and fastq_1 and fastq_2:  ## Paired-end short reads
            sample_info = ["0", fastq_1, fastq_2]
        elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
            sample_info = ["1", fastq_1, fastq_2]
        else:
            print_error("Invalid combination of columns provided!", "Line", line)

        ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2 ] }
        if sample not in sample_mapping_dict:
            sample_mapping_dict[sample] = [sample_info]
        else:
            if sample_info in sample_mapping_dict[sample]:
                print_error("Samplesheet contains duplicate rows!", "Line", line)
            else:
                sample_mapping_dict[sample].append(sample_info)

## Write validated samplesheet with appropriate columns
if len(sample_mapping_dict) > 0:
    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)
    with open(file_out, "w") as fout:
        fout.write(",".join(["sample", "single_end", "fastq_1", "fastq_2"]) + "\n")
        for sample in sorted(sample_mapping_dict.keys()):

            ## Check that multiple runs of the same sample are of the same datatype
            if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))
            ## VAOFFORD: Removed _T1 suffix to sample name
            for idx, val in enumerate(sample_mapping_dict[sample]):
                fout.write(",".join(["{}".format(sample, idx + 1)] + val) + "\n")
else:
    print_error("No entries to process!", "Samplesheet: {}".format(file_in))


"""
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
"""
