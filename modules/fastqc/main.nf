process fastqc {
    beforeScript "mkdir -p ${params.out_dir}/bcl2fastq/fastqc_out"
    publishDir "$params.out_dir/bcl2fastq/fastqc_out", mode: 'move'
    
    input:
        path bcl2

    """
    fastqc \
        $bcl2 \
        -o $params.out_dir/bcl2fastq/fastqc_out \
        -t 4
    """
}