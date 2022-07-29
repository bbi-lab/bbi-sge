#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
nextflow.enable.dsl=2


def version() {
    if (!nextflow.version.matches('>=20.04.0')) {
        println "This workflow requires Nextflow version 20.04.0 or greater -- You are running version $nextflow.version"
        exit 1
    }
}

// check whether user provides the *.config file with directory
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

process sample_sheet {

    input:
        path 

    output:
        val "*_*+*_*", emit: cigar
        val "*+*", emit: amplicon
        val "*,*", emit: experiment

    """
    python3
        ${projectDir}/bin/checking_sample_sheet.py \
        --path $params.sample_sheet \

    """
}
