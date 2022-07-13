#!/usr/bin/env nextflow

// enable dsl 2 of nextflow 
nextflow.enable.dsl=2

// need to modulate like QUANTS
// loading main.nf - *.nf - nf modules
// or loading main.nf - *.nf with in-script
// include { bcl2fastq; fastqc } from './demux.nf'

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

/*
params.sheet = 'path/test.csv'
Channel
    .fromPath(params.sheet)
    .splitCSV(header:true)
    .map{row->tuple(raw.sampleId, file(row.read1)), ....}
    .set{csv_ch}
process sample_sheet_check {
    input:
        set sample_ch 
    output:
        val "*_*+*_*", emit: cigar
        val "*+*", emit: amplicon
        val "*,*", emit: experiment
        val "", emit: prefix

    """
    python3
        ${projectDir}/bin/checking_sample_sheet.py \
        --path $params.sample_sheet \
    """
}
*/


//changing some bcl2fastq parameters to value channel
seq_dir = file( params.seq_dir )
sample_sheet = file( params.sample_sheet )

//change illumina raw data to fastq 
process bcl2fastq {
    publishDir "${params.out_dir}/bcl2fastq", mode: 'copy'
    afterScript "zgrep -c @${seq_type} *R1_* | tee read_counts.txt"

    input:
        path seq_dir
        path sample_sheet
        
    output:
        path "*_001.fastq.gz", emit: bcl2
        path "*_R1_001.fastq.gz", emit: R1
        path "*_R2_001.fastq.gz", emit: R2
        path "*"

    shell:
    '''
    bcl2fastq \
        -R !{seq_dir} -o . --sample-sheet !{sample_sheet} --no-lane-splitting --use-bases-mask Y*,I*,I*,Y* --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --interop-dir .
    rm ./Undetermined*
    '''
}


/*
** quality check of fastq
 */
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
        -f !{R1} -r !{R2} -1 ./${sample_name}.R1.fastq.gz -2 ./${sample_name}.R2.fastq.gz -A GGTTTGGAGCGAGATTGATAAAGT -B CTGAGCTCTCTCACAGCCATTTAG -M 0.1 -s ./${sample_name}.merge.fastq.gz -m 0.001 -q 20 -o 20    
    done
    '''
}

//!{params.seqprep}/${sample_name}*.fastq 
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
        python \
            !{projectDir}/bin/trimming_test.py \
            --path ${sequence} \
            --output ${sample_name}.fastq
    done
    '''
}

/*
** cDNA to gDNA
** add if statement in workflow to pass different input; output
** should work in that way
 */
/*
process rnaseq {
    input:
        path trimming

    output:
        path rna

    // only run it when rna is exsited 
    when: 
        $params.rna == "not_exist"

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
                --cDNA !{params.cDNA} \
                --amp !{params.amplicon}
                --output . \        
        else
            :
        fi
    done
    '''
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
        path trim

    output:
        path "*.sam", emit: sam
    
    //when:
    //    ${params.enrich} == "not_enrich"

    shell:
    '''
    for sequence in !{trim}; do
        s="${sequence%%.*}"; index_dot=`echo "$((${#s}))"`
        index_dash=`echo $sequence | grep -aob '-' | grep -oE '[0-9]+' | head -1`
        index_r=`echo $sequence | grep -aob 'r' | grep -oE '[0-9]+' | head -1`
        sample_name=${sequence:0:${index_dot}}
        if [[ "$sample_name" == *"r"* ]]; then
            sample_amplicon=${sample_name:0:${index_r}}
        else 
            sample_amplicon=${sequence:0:${index_dash}}
        fi
        sample_ref="/net/bbi/vol1/home/jongs2/SGE/fasta/PALB2X2.fa"
        sequence_ref="!{params.out_dir}/seqprep/trimming/$sequence"
        needleall \
            -asequence $sample_ref \
            -bsequence $sequence \
            -gapopen 10 -gapextend 0.5 \
            -outfile ./${sample_name}.sam \
            -aformat sam 
    done
    '''
}

process cigar {
    beforeScript "mkdir -p ${params.out_dir}/sam/cigar_analyzer_test"
    publishDir "${params.out_dir}/sam/cigar_analyzer_test", mode: 'move'
    
    input:
        path sam

    output:
        path "*.txt"

    //when:
    //    ${params.enrich} == "not_enrich"

    """
    python \
        $projectDir/bin/cigar_analyzer.py \
        --path ${sam} \
        --cigar ${params.cigar} \
        --output .
    """
}

process sam_to_edits {
    beforeScript "mkdir -p ${params.out_dir}/sam/edits"
    publishDir "${params.out_dir}/sam/edits", mode: 'copy'

    input:
        path sam
        
    output:
        path "*.txt", emit: edits
    
    //when:
    //    ${params.enrich} == "not_enrich"

    script:
    """
    python \
        $projectDir/bin/sam_to_edits.py \
        --amp ${params.amplicon_list} \
        --exp ${params.experiment_group} \
        --ref ${params.reference} \
        --edit ${params.edit} \
        --path ${sam}\
        --output .
    """
}

/*
** annotated variants 
*/
process annotated_variants {
    beforeScript "mkdir -p ${params.out_dir}/Final"
    publishDir "${params.out_dir}/Final", mode: 'copy'

    input:
	    path edits

    output:
        path "*.txt", emit: variants
    
    //when:
    //    ${params.enrich} == "not_enrich"

    script:
    """
    python \
        $projectDir/bin/annotated_variants.py \
        --amp ${params.amplicon_list} \
        --exp ${params.experiment_group} \
        --ref ${params.reference} \
        --edit ${params.edit} \
        --cadd ${params.cadd} \
        --clinvar ${params.clinvar} \
        --path ${edits} \
        --output .
    """
}

//take:
workflow {
    
    main:
        bcl2fastq(seq_dir, sample_sheet) // bcl2fastq - illumina to fastq
        fastqc(bcl2fastq.out.bcl2) // quality check fastqc running parrel with seqprep
        seqprep(bcl2fastq.out.bcl2, bcl2fastq.out.R1, bcl2fastq.out.R2) //
        trimming(seqprep.out.merge)
        needle(trimming.out.trim)
        //if (${params.rna} == "not_exist") {
        //    needle(trimmning.out.trim)
        //}   else {
        //    rnaseq(trimming.out.trim)
        //    needle(rnaseq.out.trim)
        //}
        cigar(needle.out.sam)
        sam_to_edits(needle.out.sam)
        annotated_variants(sam_to_edits.out.edits)
}

workflow.onComplete {
    println "bbi-sge pipeline completed at: $workflow.complete"
    println "$workflow.duration"
    println "Execution status: ${ workflow.success ? 'successed' : 'failed' }"
}

workflow.onError {
    println "Something when wrong, look up the log"
}
