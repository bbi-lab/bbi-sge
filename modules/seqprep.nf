#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
nextflow.enable.dsl=2

//changing some bcl2fastq parameters to value channel
seq_dir = file( params.seq_dir )
sample_sheet = file( params.sample_sheet )

process bcl2fastq {
    publishDir "${params.out_dir}", mode: 'copy'

    input:
        path seq_dir
        path sample_sheet
        
    output:
        path "bcl2fastq/*.fastq.gz", emit: bcl2
        path "bcl2fastq/*_R1_001.fastq.gz", emit: R1
        path "bcl2fastq/*_R2_001.fastq.gz", emit: R2

    """
    bcl2fastq \
        -R ${seq_dir} -o ./bcl2fastq  \
        --sample-sheet ${sample_sheet} --interop-dir ./bcl2fastq \
        --no-lane-splitting --use-bases-mask Y*,I*,I*,Y* \
        --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0
    """
}

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

    output:
        path "*.merge.fastq.gz", emit: merge
        path "*.R1.fastq.gz"
        path "*.R2.fastq.gz"
        path "*.txt"
    
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

process trimming {
    beforeScript "mkdir -p ${params.out_dir}/trimmning"
    publishDir "${params.out_dir}/trimming", mode: 'copy'

    input:
        path merge
    
    output: 
        path "*.fastq", emit: trim

    shell:
    """
    for sequence in !{merge}; do
        index=`echo $sequence | grep -aob '.fastq' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequence:0:${index}}
        echo $sample_name
        python \
        !{projectDir}/bin/trimming.py \
        --path !{params.out_dir}/seqprep/merge/!{merge} \
        --output ./${sample_name}.fastq
    done
    """
}



process needle {
    beforeScript "mkdir -p ${params.out_dir}/sam"
    publishDir "$params.out_dir/sam", mode: 'copy'

    input:
        path trimming
    //if (params.rna == true) {
    //    path trimming_and_rnaseq
    //}   else {
    //    path trimming
    //}
    // going to pass seqeuences as inputs (.fasta)

    //working_dir = file("/net/bbi/vol1/home/jongs2/needle_test")
    //File[] listOfFiles = working_dir.listFiles()

    //index_dot=`echo $sequence | grep -aob '\.' | grep -oE '[0-9]+' | head -1`
    //maybe you need to escape # with \
    output:
        path "*.sam", emit: sam
    
    shell:
    '''
    for sequence in !{merge}; do
        s="${sequence%%.*}"; index_dot=`echo "$((${#s}))"`
        index_dash=`echo $sequence | grep -aob '-' | grep -oE '[0-9]+' | head -1`
        index_r=`echo $sequence | grep -aob 'r' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequence:0:${index_dot}}
        if [[ "$sample_name" == *"r"* ]]; then
            sample_amplicon=${sample_name:0:${index_r}}
        else 
            sample_amplicon=${sequence:0:${index_dash}}
        fi
        sample_ref="/net/bbi/vol1/home/jongs2/SGE/fasta/${sample_amplicon}.fa"
        sequence_ref="!{params.out_dir}/seqprep/merge/!{merge}"
        needleall \
            -asequence $sample_ref \
            -bsequence $sequence_ref \
            -gapopen 10 -gapextend 0.5 \
            -outfile ./${sample_name}.sam \
            -aformat sam 
    done
    '''
}

workflow {
    main:
        bcl2fastq(seq_dir, sample_sheet)
        seqprep(bcl2fastq.out.R1, bcl2fastq.out.R2)
        trimmning(seqprep.out.merge)
        needle(trimming.out.trim)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}