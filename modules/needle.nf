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

params.testing_out = "/net/bbi/vol1/home/jongs2/test/20220110_Nextseqrun_PracticeRunthroughX5DX11/nobackup/fastq/Seqprep/merged/no_Ns"
process needle {
    beforeScript "mkdir -p ${params.out_dir}/sam"
    publishDir "$params.out_dir/sam", mode: 'copy'

    //input:
    //if (params.rna == true) {
    //    path trimming_and_rnaseq
    //}   else {
    //    path trimming
    //}
    // going to pass seqeuences as inputs (.fasta)

    //working_dir = file("/net/bbi/vol1/home/jongs2/needle_test")
    //File[] listOfFiles = working_dir.listFiles()

    //maybe you need to escape # with \
    output:
        path "*.sam", emit: sam
    
    shell:
    '''
    for sequence in "!{params.testing_out}/*.fastq"; do
        echo $sequence
        index_dot=`echo $sequence | grep -aob '\.' | grep -oE '[0-9]+' | head -1`
        index_dash=`echo $sequence | grep -aob '-' | grep -oE '[0-9]+' | head -1`
        index_r=`echo $sequence | grep -aob 'r' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequence:0:${index_dot}}
        if [[ "$sample_name" == *"r"* ]]; then
            sample_amplicon=${sample_name:0:${index_r}}
        else 
            sample_amplicon=${sequence:0:${index_dash}}
        fi
        sample_ref="/net/bbi/vol1/home/jongs2/SGE/fasta/${sample_amplicon}.fa"
        sequence_ref="!{params.out_dir}/seqprep/merge/$sequence"
        needleall \
            -asequence $sample_ref \
            -bsequence $sequence_ref \
            -gapopen 10 -gapextend 0.5 \
            -outfile ./${sample_name}.sam \
            -aformat sam 
    done
    '''
}

workflow {
    main:
        needle()
}

workflow.onComplete {
    println "bbi-sge pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}
