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

// some defined parameters
params.help = false

// help message
if (params.help) {
    log.info
    '''
    -----------------------------------------------------------------------------------------------
    Required parameters
        params.out_dir             - output file dir
        params.seq_dir             - sequencing reading dir
        params.sample_sheet        - sample sheet dir
        params.seq_type            - "NS" for Nextseq / "M" for Miseq
        params.rna                 - check box for rna existed or not (True = existed; False = not existed)
        params.ref_dir             -
        params.cigar_comparison    -
        params.amplincon_list      -
        params.experiment_grouping -
        params.prefix              -
    -----------------------------------------------------------------------------------------------
    Optional parameters
        Unless you change the parameters in specific in here, it will take the parameters that were setup in main.nf
        -------------------------------------------------------------------------------------------
        seqprep
            
        -------------------------------------------------------------------------------------------
        needle
    -----------------------------------------------------------------------------------------------
    '''
    exit 1
}

//changing some bcl2fastq parameters to value channel
seq_dir = file( params.seq_dir )
sample_sheet = file( params.sample_sheet )

/*
** change illumina data to fastq 
 */
process bcl2fastq {
    // replaceFisrt, replaceAll
    // wordStartsWithS = ~/(?i)\s+Gr\w+/
    // delete Ss
    publishDir "${params.out_dir}", mode: 'copy'//, saveAs: { *.fastq.gz -> "*."}
    // afterScript vs. innerscript (not working) vs. def. vs. process
    afterScript "zgrep -c @${seq_type} bcl2fastq/*R1_* | tee read_counts.txt"
    // channel filter??
    //myFile = file('Undetermined*')
    //result = myFile.delete()
    //afterScript "rm Undetermined*"

    input:
        path seq_dir
        path sample_sheet
        
    output:
        path "bcl2fastq/*.fastq.gz", emit: bcl2
        path "bcl2fastq/*_R1_001.fastq.gz", emit: R1
        path "bcl2fastq/*_R2_001.fastq.gz", emit: R2
        //[pre, post, lib, neg]
        //path (to pass non-*.fastq.gz)
        //path "bcl2fastq"

    """
    bcl2fastq \
        -R ${seq_dir} -o ./bcl2fastq  \
        --sample-sheet ${sample_sheet} --interop-dir ./bcl2fastq \
        --no-lane-splitting --use-bases-mask Y*,I*,I*,Y* \
        --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0
    """
}

/*
** quality check of fastq
 */

process fastqc {
    //beforeScript vs. publishDir vs. both
    beforeScript "mkdir -p ${params.out_dir}/fastqc_out"
    publishDir "${params.out_dir}/fastqc_out", mode: 'move'
    cache 'lenient'
    
    input:
        path bcl2

    """
    fastqc \
        ${bcl2} -o ${params.out_dir}/fastqc_out -t 8
    """
}

/*
process seqprep {
    beforeScript "mkdir -p ${params.out_dir}/seqprep"
    beforeScript "mkdir -p ${params.out_dir}/seqprep/R1"
    beforeScript "mkdir -p ${params.out_dir}/seqprep/R2"
    beforeScript "mkdir -p ${params.out_dir}/seqprep/merge"
    publishDir "${params.out_dir}/seqprep/R1", mode: 'copy'//, saveAs: "*.R1.fastq.gz"
    publishDir "${params.out_dir}/seqprep/R2", mode: 'copy'//, saveAs: "*.R2.fastq.gz"
    publishDir "${params.out_dir}/seqprep/merge", mode: 'copy'//, saveAs: "*.merged.fastq.gz"
    afterScript "zgrep -c @${seq_type} *.fastq.gz >> seqprep_read_counts.txt"

    input:
        path R1
        path R2

    output:
        path "seqprep/merge/*.merge.fastq.gz", emit: merge
        path "seqprep/merge/*.R1.fastq.gz"
        path "seqprep/merge/*.R2.fastq.gz"
        
    
    script:
    """
    for i in "${R1}"
        do
            index_of_first_under = i.find('_')
            sample_name = i[:index_of_first_under]
            $projectDir/SeqPrep/./SeqPrep \
                -f "${R1}" -r "${R2}"\
                -1 /seqprep/R1/"${sample_name}".R1.fastq.gz
                -2 /seqprep/R2/"${sample_name}".R2.fastq.gz \
                -s /seqprep/merge/"${sample_name}".merge.fastq.gz \
                -A GGTTTGGAGCGAGATTGATAAAGT -B CTGAGCTCTCTCACAGCCATTTAG \
                -M 0.1 -m 0.001 -q 20 -o 20)
        done
    """
}
*/
/*
process trimming {
    beforeScript "mkdir -p ${params.out_dir}/trimmning"
    publishDir "${params.out_dir}/trimming", mode: 'copy'

    input:
        path merge
    
    output: 
        path trimming

    """
	python \
		$projectDir/bin/trimming.py --path {params.out_dir}
    """
}
*/

workflow {
    main:
        bcl2fastq(seq_dir, sample_sheet)
        fastqc(bcl2fastq.out.bcl2)
        //seqprep(bcl2fastq.out.R1, bcl2fastq.out.R2)
        //trimming(seqperp.out.merge)
}

workflow.onComplete {
    println "bbi-sge pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}
