#!/usr/bin/env nextflow

include { sam_to_edits } from '../modules/edits.nf'
workflow EDITS {
    main:
        sam_to_edits()
}

