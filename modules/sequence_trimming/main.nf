process trimming {
    beforeScript "mkdir -p ${params.out_dir}/trimming"
    publishDir "${params.out_dir}/trimming", mode: 'copy'

    input:
        path merge
    
    output: 
        path "*.fastq", emit: trimming

    shell:
    myDir = file("/net/bbi/vol1/home/jongs2/test/20220110_Nextseqrun_PracticeRunthroughX5DX11/nobackup/fastq/Seqprep/merged")
    printnln(myDir)
    listOfFiles = myDir.list()
    '''
    for sequence in !{listOfFiles}; do
        echo $sequence
        index=`echo $sequence | grep -aob '.fastq' | grep -oE '[0-9]+' | head -1`
        echo $index
        sample_name=${sequence:0:${index}}
        python \
            !{projectDir}/bin/trimming.py \
            --path !{params.testing} \
            --output ./${sample_name}.fastq
    done
    '''
}