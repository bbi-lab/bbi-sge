#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2

include { cigar } from '../modules/cigar.nf'
workflow CIGAR{
    take:
        reads // file ( .sam files after alignment )
        amp // file ( .fastq file for amplicon sequence )
        exp // experiment list
        ref // file ( reference directory for amplicon )
        edit // file ( directory of edit.txt )
        output // file ( output directory )


    main:
        cigar(reads, )
        cigar

    emit:
        
}

include { sam_to_edits } from '../modules/edits.nf'
workflow EDITS {
    take:
        reads // file ( .sam files after alignment )

    main:
        sam_to_edits()
}

include { annotated_variants } from '../modules/annotation.nf'
workflow ANNOTATION {
    take:
        edits // file ( .txt files after edits ) ()

    main:
        annotated_variants()
}

