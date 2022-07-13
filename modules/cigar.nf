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


process cigar {
    
    beforeScript "mkdir -p ${params.out_dir}/sam/cigar_analyzer"
    publishDir "${params.out_dir}/sam/cigar_analyzer", mode: 'copy'
    
    //input:
    //just need to put input (*.sam) and replace {params.testing_sam}

    output:
        path "*.txt"
    """
    python \
        $projectDir/bin/cigar_analyzer.py \
        --path ${params.testing_sam} \
        --cigar ${params.cigar} \
        --output . 
    head -15 cigar_counts/*.txt  >> cigar_counts/combined_cigar_counts.txt
    """
}

workflow {
    main:
        cigar()
}

workflow.onComplete {
    println "bbi-sge pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}
