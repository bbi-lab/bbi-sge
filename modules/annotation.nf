// check whether user provides the *.config file with directory
//if (!params.seq_dir || !params.out_dir || !params.sample_sheet) {
//    exit 1, 
//	"Must include config file using -c <config name>.config that includes seq_dir, out_dir, and sample_sheet."
//}

//changing some bcl2fastq parameters to value channel
        //--amp ${params.amplicon_list}
        //--exp ${params.experiment_group}
params.testing_variants = "/net/bbi/vol1/home/jongs2/test/20220110_Nextseqrun_PracticeRunthroughX5DX11/nobackup/fastq/Seqprep/merged/no_Ns/sam/variant_counts_no_thresh"
params.testing_variants2 = "/mnt/c/Users/shinj/Desktop/CAVA/CAVA_Production/Pipeline/20220110_Nextseqrun_PracticeRunthroughX5DX11/nobackup/fastq/Seqprep/merged/no_Ns/sam/variant_counts_no_thresh"
process annotated_variants {
    beforeScript "mkdir -p ${params.out_dir}/edits/annotated"
    publishDir "${params.out_dir}/edits/annotated", mode: 'copy'

    output:
        path "*.txt", emit: annotated
    
    script:
    """
    python \
        $projectDir/bin/annotated_variants_test.py \
        --amp ${params.amplicon_list} \
        --exp ${params.experiment_group} \
        --path ${params.testing_variants2} \
        --output .
    """
}