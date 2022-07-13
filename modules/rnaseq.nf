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

process rnaseq {
    //beforeScript "mkdir -p ${params.out_dir}/sam"
    publishDir "${params.out_dir}/test", mode: 'copy' 
    //input:
    //    path cDNA

    //output:
    //    path gDNA

    // only run it when rna is existed
    when: 
        params.rna == 'exist'



    shell:
    '''
    for i in "!{params.out_dir}/seqprep/merge/*.fastq"; do
        if [[ ${i} == *.fastq ] && [ 'rna' in ${i}]]; then
            print i
            index_dot=`echo $sequence | grep -aob '.' | grep -oE '[0-9]+' | head -1`
            index_r=`echo $sequence | grep -aob 'r' | grep -oE '[0-9]+' | head -1`
            index_extension=`echo $sequence | grep -aob '.fastq' | grep -oE '[0-9]+' | head -1`
            amp=${i:0:${index_r}}; sample_name=${i:0:${index_dot}}; before_extension_name=${i:0:${index_extension}}
            python \
                !{projectDir}/bin/rnaseq.py \
                --path !{params.testing_sam} \
                --output . \
                --amp $amp
        else
            pass
        fi
    done
    '''
}

workflow {
    main:
        rnaseq()
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}