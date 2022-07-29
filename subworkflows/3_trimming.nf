#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2


include { trimming  } from '../modules/sequence_trimming/main.nf'
workflow SEQUENCE_TRIM{
    main:
        trimming()
}

include { rna_trimming } from '../modules/rna_trimming/main.nf'
workflow RNA_TRIM {
    main:
        rna_trimming()
}


