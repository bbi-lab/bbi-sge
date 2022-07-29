#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2


}

include { needleall } from '../modules/needleall/main.nf'
/*
//
*/
workflow ALIGNMENT {
    main:
        needleall()
}

