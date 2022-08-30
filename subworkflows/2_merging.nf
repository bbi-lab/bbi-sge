#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2

// loading seqprep module 
include { seqprep } from '../modules/seqprep/main.nf'
/*
** ... merging and also a bit of adapter trimming
*/
workflow MERGING {
    take: 
        directory

    main: 
        seqprep(R1, R2)
        ch_merge = seqprep.out.merge

    emit:
    reads = merge
}