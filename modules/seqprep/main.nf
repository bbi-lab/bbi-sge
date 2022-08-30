#!/usr/bin/env nextflow

process seqprep {
    beforeScript "mkdir -p ${params.out_dir}/seqprep/R1"
    beforeScript "mkdir -p ${params.out_dir}/seqprep/R2"
    beforeScript "mkdir -p ${params.out_dir}/seqprep/merge"
    publishDir "${params.out_dir}/seqprep/R1", mode: 'copy', pattern: "*.R1.fastq.gz"
    publishDir "${params.out_dir}/seqprep/R2", mode: 'copy', pattern: "*.R2.fastq.gz"
    publishDir "${params.out_dir}/seqprep/merge", mode: 'copy', pattern: "*.merge.fastq.gz"
    afterScript "zgrep -c @${seq_type} *.fastq.gz >> seqprep_read_counts.txt"

    input:
        path directory
        //val prefix

    output:
        path "*.merge.fastq.gz", emit: merge
        path "*.R1.fastq.gz"
        path "*.R2.fastq.gz"
        path "*.txt"

    script:
        def foo = R1.collect{"-I $it"}.join(' ')
        def bar = R2.collect{"-I $it"}.join(' ')
    
    shell:
    '''
    for sequences in !{directory}; do
        !{projectDir}/SeqPrep/./SeqPrep \
            -f !{directory}/*_R1_001.fastq.gz \
            -r !{directory}/*_R2_001.fastq.gz \
            -1 ./test.R1.fastq.gz \
            -2 ./test.R2.fastq.gz \
            -s ./test.merge.fastq.gz \
            -A GGTTTGGAGCGAGATTGATAAAGT \
            -B CTGAGCTCTCTCACAGCCATTTAG \
            -M 0.1 -m 0.001 -q 20 -o 20        
    done
    '''
}