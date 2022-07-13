#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
nextflow.enable.dsl=2

//changing some bcl2fastq parameters to value channel
seq_dir = file( params.seq_dir )
sample_sheet = file( params.sample_sheet )

/*
** change illumina raw data to fastq 
 */
process bcl2fastq {

    publishDir "${params.out_dir}/bcl2fastq", mode: 'copy'
    afterScript "zgrep -c @${seq_type} *R1_* | tee read_counts.txt"

    input:
        path seq_dir
        path sample_sheet
        
    output:
        path "${params.prefix}*.fastq.gz", emit: bcl2
        path "${params.prefix}*_R1_001.fastq.gz", emit: R1
        path "${params.prefix}*_R2_001.fastq.gz", emit: R2
        path "*"

    """
    bcl2fastq \
        -R ${seq_dir} -o . \
        --sample-sheet ${sample_sheet} --interop-dir . \
        --no-lane-splitting --use-bases-mask Y*,I*,I*,Y* \
        --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0
    """
}

process fastqc {
    beforeScript "mkdir -p ${params.out_dir}/bcl2fastq/fastqc_out"
    publishDir "$params.out_dir/bcl2fastq/fastqc_out", mode: 'move'
    
    input:
        path bcl2

    """
    fastqc \
        $bcl2 -o $params.out_dir/bcl2fastq/fastqc_out -t 4
    """
}
workflow {
    main:
        bcl2fastq()
        fastqc(bcl2fastq.out.bcl2)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}