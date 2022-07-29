process rna_trimming {
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