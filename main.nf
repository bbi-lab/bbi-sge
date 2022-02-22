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

/*
process sample_sheet_check {
    publishDir path: "${params.out_dir}", pattern: "SampleSheet.csv", mode: 'copy'
    publishDir path: "${params.out_dir}", pattern: "SampleMap.csv", mode: 'copy'
    publishDir path: "${params.out_dir}", pattern: "GarnettSheet.csv", mode: 'copy'
    publishDir path: "${params.out_dir}/sample_id_maps", pattern: "*_SampleIDMap.csv", mode: 'copy'

    input:
        file insamp from Channel.fromPath(params.generate_samplesheets)

    output:
        file "*Sheet.csv"
        file "SampleMap.csv" optional true
        file "*_SampleIDMap.csv" optional true

    when:
        params.generate_samplesheets != 'no_input'


    """
    generate_sample_sheets.py $params.generate_samplesheets
    """
}
*/

//changing some bcl2fastq parameters to value channel
seq_dir = file( params.seq_dir )
sample_sheet = file( params.sample_sheet )

/*
** change illumina raw data to fastq 
 */
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

/*
** quality check of fastq
 */
/*
process fastqc {
    beforeScript "mkdir -p ${params.out_dir}/bcl2fastq/fastqc_out"
    publishDir "$params.out_dir/bcl2fastq/fastqc_out", mode: 'move'
    
    input:
        path bcl2

    """
    fastqc \
        $bcl2 -o $params.out_dir/bcl2fastq/fastqc_out -t 4
    """
}
*/

/*
** https://github.com/jstjohn/SeqPrep
** SeqPrep is a program to merge paired end Illumina reads that are overlapping into a single longer read.
 */
process seqprep {
    beforeScript "mkdir -p ${params.out_dir}/seqprep/R1"
    beforeScript "mkdir -p ${params.out_dir}/seqprep/R2"
    beforeScript "mkdir -p ${params.out_dir}/seqprep/merge"
    publishDir "${params.out_dir}/seqprep/R1", mode: 'copy', pattern: "*.R1.fastq.gz"
    publishDir "${params.out_dir}/seqprep/R2", mode: 'copy', pattern: "*.R2.fastq.gz"
    publishDir "${params.out_dir}/seqprep/merge", mode: 'copy', pattern: "*.merge.fastq.gz"
    afterScript "zgrep -c @${seq_type} *.fastq.gz >> seqprep_read_counts.txt"

    input:
        path R1
        path R2

    output:
        path "*.merge.fastq.gz", emit: merge
        path "*.R1.fastq.gz"
        path "*.R2.fastq.gz"
        path "*.txt"
    
    shell:
    '''
    for sequences in !{R1}; do
        index=`echo $sequences | grep -aob '_' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequences:0:${index}}
        !{projectDir}/SeqPrep/./SeqPrep \
            -f !{params.out_dir}/bcl2fastq/!{R1} \
            -r !{params.out_dir}/bcl2fastq/!{R2} \
            -1 ./${sample_name}.R1.fastq.gz \
            -2 ./${sample_name}.R2.fastq.gz \
            -s ./${sample_name}.merge.fastq.gz \
            -A GGTTTGGAGCGAGATTGATAAAGT -B CTGAGCTCTCTCACAGCCATTTAG \
            -M 0.1 -m 0.001 -q 20 -o 20        
    done
    '''
}

process trimming {
    beforeScript "mkdir -p ${params.out_dir}/seqprep/trimming"
    publishDir "${params.out_dir}/seqprep/trimming", mode: 'copy'

    input:
        path merge
    
    output: 
        path "*.fastq", emit: trim

    shell:
    '''
    for sequence in !{merge}; do
        index=`echo $sequence | grep -aob '.fastq' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequence:0:${index}}
        echo $sample_name
        python \
            !{projectDir}/bin/trimming.py \
            --path !{params.out_dir}/seqprep/merge/$sequence \
            --output ${sample_name}.fastq
    done
    '''
}

/*
** cDNA to gDNA
 */
/*
process rnaseq {
    input:
        path trimming

    output:
    

    // only run it when rna is exsited 
    when: 
        $params.rna != "not_exist"

    """
    python \
    """
}
*/

/*
** alignment to sam
** need to figure out to bam
 */
process needle {
    beforeScript "mkdir -p ${params.out_dir}/sam"
    publishDir "$params.out_dir/sam", mode: 'copy'

    input:
        path trimming

    output:
        path "*.sam", emit: sam
    
    shell:
    '''
    for sequence in !{trimming}; do
        s="${sequence%%.*}"; index_dot=`echo "$((${#s}))"`
        index_dash=`echo $sequence | grep -aob '-' | grep -oE '[0-9]+' | head -1`
        index_r=`echo $sequence | grep -aob 'r' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequence:0:${index_dot}}
        if [[ "$sample_name" == *"r"* ]]; then
            sample_amplicon=${sample_name:0:${index_r}}
        else 
            sample_amplicon=${sequence:0:${index_dash}}
        fi
        sample_ref="/net/bbi/vol1/home/jongs2/SGE/fasta/${sample_amplicon}.fa"
        sequence_ref="!{params.out_dir}/seqprep/trimming/$sequence"
        needleall \
            -asequence $sample_ref \
            -bsequence $sequence_ref \
            -gapopen 10 -gapextend 0.5 \
            -outfile ./${sample_name}.sam \
            -aformat sam 
    done
    '''
}

process cigar {
    
    beforeScript "mkdir -p ${params.out_dir}/sam/cigar_analyzer"
    publishDir "${params.out_dir}/sam/cigar_analyzer", mode: 'move'
    
    input:
        path sam
    //just need to put input (*.sam) and replace {params.testing_sam}

    output:
        path "*.txt"

    when:
        ${params.enrich} != "not_enrich"

    """
    python \
        $projectDir/bin/cigar_analyzer.py \
        --path ${params.testing_sam} \
        --cigar ${params.cigar} \
        --output . 
    head -15 cigar_counts/*.txt  >> cigar_counts/combined_cigar_counts.txt
    """
}
/*
process sam_to_edits {
    beforeScript "mkdir -p ${params.out_dir}/edits"
    publishDir "${params.out_dir}/edits", mode: 'copy'

    //input:
        //path sam
        
    output:
        path "*.txt", emit: edits
    
    when:
        ${params.enrich} != "not_enrich"

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
*/

/*
** annotated variants 
*/
/*
process annotated_variants {
    beforeScript "mkdir -p ${params.out_dir}/annotated"
    publishDir "${params.out_dir}/annotated", mode: 'copy'

    output:
        path "*.txt", emit: annotated
    
    when:
        ${params.enrich} != "not_enrich"

    script:
    """
    python \
        $projectDir/bin/annotated_variants.py \
        --amp ${params.amplicon_list} \
        --exp ${params.experiment_group} \
        --path ${params.testing_sam} \
        --output .
    """
}
*/

workflow {
    //take: 
    //  fastqc = Channel.fromPath(bcl2fastq.out_bcl2)

    main:
        bcl2fastq(seq_dir, sample_sheet) // bcl2fastq - illumina to fastq
        //fasqc(bcl2fastq.out.bcl2) // quality check fastqc running parrel with seqprep
        seqprep(bcl2fastq.out.R1, bcl2fastq.out.R2) //
        trimming(seqprep.out.merge)
        needle(trimming.out.trim)
        //if (${params.rna} == "exist") {
        //    rnaseq(trimming.out.trim)
        //    needle(rnaseq.out.gDNA)
        //}   else {
        //    needle(trimmning.out.trim)
        //}
        //if (${params.enrich} == "enrich")
        cigar(needle.out.sam)
        //sam_to_edits(needle.out.sam)
        //annotated_variants(sam_to_edits.out)
}

workflow.onComplete {
    println "bbi-sge pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'successed' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}
