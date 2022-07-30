process bcl2fastq {
    publishDir "${params.out_dir}/bcl2fastq", mode: 'copy'
    afterScript "zgrep -c @${seq_type} *R1_* | tee read_counts.txt"

    input:
        path seq_dir
        path sample_sheet
        
    output:
        path "${params.prefix}*.fastq.gz", emit: bcl2
        path "${params.prefix}*_R1_001.fastq.gz", emit: R1
        path "${params.prefix}*_R2_001.fastq.gz", emit: R2
        path "*"

    """
    bcl2fastq \
        -R ${seq_dir} \
        -o . \
        --sample-sheet ${sample_sheet} \
        --interop-dir . \
        --no-lane-splitting \
        --use-bases-mask Y*,I*,I*,Y* \
        --minimum-trimmed-read-length 0 \
        --mask-short-adapter-reads 0
    """
}