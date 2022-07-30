#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2


include { annotated_variants } from '../modules/annotation.nf'
workflow ANNOTATION {


    main:
        annotated_variants()
}

