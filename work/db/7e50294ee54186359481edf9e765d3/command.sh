#!/bin/bash -ue
bcl2fastq         -R 211123_NB552332_0205_AHCGCFAFX3 -o ./bcl2fastq          --sample-sheet SampleSheet.csv --interop-dir ./bcl2fastq         --no-lane-splitting --use-bases-mask Y*,I*,I*,Y*         --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0
