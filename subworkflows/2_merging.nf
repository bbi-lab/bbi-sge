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
        R1     // file ( r1 for demux )
        R2     // file ( r2 for demux)
        prefix // val  ( stored prefix from sampe_sheet )

    main: 
        seqprep(R1, R2)
        ch_merge = seqprep.out.merge

    //emit:
    reads = ch_merge
}