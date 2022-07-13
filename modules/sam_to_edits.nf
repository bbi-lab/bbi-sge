#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
nextflow.enable.dsl=2

// check the version of Nextflow
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

//changing some bcl2fastq parameters to value channel
        //--amp ${params.amplicon_list}
        //--exp ${params.experiment_group}

// for right now do not set up input and setup path as files from greg
params.testing_sam_to_edits = "/net/bbi/vol1/home/jongs2/test/20220110_Nextseqrun_PracticeRunthroughX5DX11/nobackup/fastq/Seqprep/merged/no_Ns/sam"
process sam_to_edits {
    beforeScript "mkdir -p ${params.out_dir}/edits"
    publishDir "${params.out_dir}/edits", mode: 'copy'

    //input:
        //path sam
        
    output:
        path "*.txt", emit: edits
    
    script:
    """
    python \
        $projectDir/bin/sam_to_edits.py \
        --amp ${params.amplicon_list} \
        --exp ${params.experiment_group} \
        --path ${params.testing_sam_to_edits} \
        --output .
    """
}

workflow {
    main:
        sam_to_edits()
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}