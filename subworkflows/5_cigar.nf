#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2

include { cigar } from '../modules/cigar.nf'
workflow CIGAR{
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
