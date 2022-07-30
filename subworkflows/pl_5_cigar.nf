#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2

include { cigar } from '../modules/cigar.nf'
workflow CIGAR{
    main:
        cigar()
}
