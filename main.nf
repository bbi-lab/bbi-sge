#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
nextflow.enable.dsl=2

// check the version of Nextflow
// if the condition is not met, exit
def version() {
    if (!nextflow.version.matches('>=20.04.0')) {
        println "This workflow requires Nextflow version 20.04.0 or greater -- You are running version $nextflow.version"
        exit 1
    }
}

// check whether user provides the *.config file with directory
// if the condition is not met, exit
if (!params.seq_dir || !params.out_dir || !params.sample_sheet) {
    exit 1, 
	"Must include config file using -c <config name>.config that includes seq_dir, out_dir, and sample_sheet."
}


// some defined parameters
params.help = false

// help message
if (params.help) {
    log.info
    '''
    -----------------------------------------------------------------------------------------------
    Required parameters
    These parameters are required to run the pipeline
        params.out_dir             - the path directory of your output files
        params.seq_dir             - the path diretory of your sequence files
        params.sample_sheet        - the path of your sample sheet (.csv)
        params.seq_type            - "NS" for Nextseq / "M" for Miseq
        params.rna                 - check box for rna existed or not (default: not_exist; change to whatever you want it)
        params.enrich              - check box for enrich step or entire step (default: not_enrich; change to whatever you want it)
        params.ref_dir             - the mother path directory of your reference sequences for cigar, amplicon, cadd, clinvar and others

        params.cigar               - 
        params.amplincon_list      -
        params.experiment_group    -
        
    -----------------------------------------------------------------------------------------------
    DO NOT NEED TO CARE ABOUT THESE FOR NOW
    -----------------------------------------------------------------------------------------------
        params.name - integrated with sample name of sample sheet or after bcl2fastq for other variable names (cigar, amplicon, etc.)
    -----------------------------------------------------------------------------------------------
    cadd 
        params.chromosome   -
        params.basepairs    - number of base pairs you are searching for (ex.50, 120, 200, 400; wheatever)
        params.starting     - staring base pair position
        params.ending       - ending base pair position
        
        *basepar is

    clinvar

    -----------------------------------------------------------------------------------------------
    '''
    exit 1
}

// seperate demux and non-demux
// seqprep change... which should not give any errors
// if {== run}
// else 
// give file directory ==
// or overwrite input



// raw data -> sam
include { SAMPLE } from './subworkflows/0_sample_sheet.nf'            // check_sample_sheet input and pass variables accordingly
include { DEMUX; DEMUX_QC } from './subworkflows/1_demux.nf'                // demux and raw_qc using 'bcl2fastq' and 'fastqc'
include { MERGING } from "./subworkflows/2_seqprep.nf"                      // merging and adapter trimming using 'seqprep'
include { SEQUENCE_TRIM; RNA_TRIM } from './subworkflows/3_trimming.nf'     // some other trimming - Ns, rna
//include { ALIGNMENT } from './subworkflows/4_alignment.nf'                  // alignment using 'needleall'

// analysis after sam
//include { CIGAR; EDITS; ANNOTATION } from './subworkflows/5_post_alignment' // 




workflow {

    // defining the critical inputs
    seq_dir = file( "${params.out_dir}/bcl2fastq" )

    


// comment out from here; not yet tested, placehodler for now
    if (params.merging == 'seqprep') {
        MERGING(seq_dir)
        merge = MERGING.out.merge
    } else {
        merge = reads
    }


}

workflow.onComplete {
    println "bbi-sge pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'successed' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}