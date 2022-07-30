#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2


// including needleall module
include { needleall } from '../modules/needleall/main.nf'

/*
** needleall alignment
** needle... link
** 
*/
workflow ALIGNMENT {
    take:
        reads

    main:
        needleall()

    emit:

}

// loading fastqc module
include { fastqc } from '../modules/fastqc/main.nf'
/*
** alignment fastqc 
*/
workflow ALIGNMENT_QC {
    take: 
        reads

    main: 
        fastqc(reads)
}

