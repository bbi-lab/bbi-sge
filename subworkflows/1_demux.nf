#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
//nextflow.enable.dsl=2

// loading bcl2fastq module
include { bcl2fastq } from '../modules/bcl2fastq/main.nf'
/*
** demux raw bcl2fastq data to fastq format
** reads - all 
** r1 -
** r2 -
** detailed in .yaml 
** link on here or .yaml?
*/
workflow DEMUX {
    take: 
        seq_dir      // file ( file directory of illumina raw sequences )
        sample_sheet // file ( .csv of sample_sheet )
        //prefix       // val  ( stored prefix from sampe_sheet )

    main: 
        bcl2fastq(seq_dir, sample_sheet)
        bcl2 = bcl2fastq.out.bcl2
        r1 = bcl2fastq.out.R1
        r2 = bcl2fastq.out.R2

    emit:
    reads = bcl2
    r1_reads = r1
    r2_reads = r2
}


// loading fastqc module
include { fastqc } from '../modules/fastqc/main.nf'
/*
** raw sequencing fastqc 
** qc...? 
** emit required??? only if I will use it...
*/
workflow DEMUX_QC {
    take: 
        reads

    main: 
        fastqc(reads)
}