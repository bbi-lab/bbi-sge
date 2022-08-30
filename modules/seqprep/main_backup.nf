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
        path R1
        path R2
        val prefix

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
    for sequences in !{R1}; do
        index=`echo $sequences | grep -aob '_' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequences:0:${index}}
        !{projectDir}/SeqPrep/./SeqPrep \
            -f !{params.out_dir}/bcl2fastq/!{R1} \
            -r !{params.out_dir}/bcl2fastq/!{R2} \
            -1 ./${sample_name}.R1.fastq.gz \
            -2 ./${sample_name}.R2.fastq.gz \
            -s ./${sample_name}.merge.fastq.gz \
            -A GGTTTGGAGCGAGATTGATAAAGT -B CTGAGCTCTCTCACAGCCATTTAG \
            -M 0.1 -m 0.001 -q 20 -o 20        
    done
    '''
}