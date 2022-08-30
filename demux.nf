#!/usr/bin/env nextflow

/*
** Seperated the demux process from all other prcoess due to two reasons
** 1. There are few times to manipulate the sequences after demux and it is easier with this way.
** 2. It is easier to pass the variables for later in pipeline.
*/


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

//changing some bcl2fastq parameters to value channel
seq_dir = file( params.seq_dir )
sample_sheet = file( params.sample_sheet )

process bcl2fastq {

    publishDir "${params.out_dir}/bcl2fastq", mode: 'copy'//, saveAs: { *.fastq.gz -> "*."}
    afterScript "zgrep -c @${seq_type} *R1_* | tee read_counts.txt"
    //afterScript "rm Undetermined*"

    input:
        path seq_dir
        path sample_sheet

    output:
        path "*.fastq.gz", emit: bcl2
        path "*_R1_001.fastq.gz", emit: R1
        path "*_R2_001.fastq.gz", emit: R2
        path "*.txt"

    """
    bcl2fastq \
        -R ${seq_dir} -o . \
        --sample-sheet ${sample_sheet} --interop-dir . \
        --no-lane-splitting --use-bases-mask Y*,I*,I*,Y* \
        --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0
    """
}

workflow {
    bcl2fastq (seq_dir, sample_sheet)
}

workflow.onComplete {
    println "bbi-sge pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'successed' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}