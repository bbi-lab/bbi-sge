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
