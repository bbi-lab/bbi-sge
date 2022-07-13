#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
nextflow.enable.dsl=2

// check the version of Nextflow
def version() {
    if (!nextflow.version.matches('>=20.04.0')) {
        println "This workflow requires Nextflow version 20.04.0 or greater -- You are running version $nextflow.version"
        exit 1
    }
}

// check whether user provides the *.config file with directory
if (!params.seq_dir || !params.out_dir || !params.sample_sheet) {
    exit 1, 
	"Must include config file using -c <config name>.config that includes seq_dir, out_dir, and sample_sheet."
}

process trimming {
    beforeScript "mkdir -p ${params.out_dir}/trimming"
    publishDir "${params.out_dir}/trimming", mode: 'copy'

    input:
        path merge
    
    output: 
        path "*.fastq", emit: trimming

    shell:
    myDir = file("/net/bbi/vol1/home/jongs2/test/20220110_Nextseqrun_PracticeRunthroughX5DX11/nobackup/fastq/Seqprep/merged")
    printnln(myDir)
    listOfFiles = myDir.list()
    '''
    for sequence in !{listOfFiles}; do
        echo $sequence
        index=`echo $sequence | grep -aob '.fastq' | grep -oE '[0-9]+' | head -1`
        echo $index
        sample_name=${sequence:0:${index}}
        python \
            !{projectDir}/bin/trimming.py \
            --path !{params.testing} \
            --output ./${sample_name}.fastq
    done
    '''
}

workflow {
    main:
        trimming()
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}