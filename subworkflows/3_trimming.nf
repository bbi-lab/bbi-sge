#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2

// maybe adapter trimming ( seqprep ), removing Ns ( sequence_trimming ) and for one...?

// including trimming module
include { trimming  } from '../modules/sequence_trimming/main.nf'

/*
** removing Ns 
*/
workflow SEQUENCE_TRIM{
    main:
        trimming()
}

// including rna_trimming module 
include { rna_trimming } from '../modules/rna_trimming/main.nf'

/*
**
*/
workflow RNA_TRIM {
    main:
        rna_trimming()
}


