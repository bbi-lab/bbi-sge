process cigar {
    beforeScript "mkdir -p ${params.out_dir}/sam/cigar_analyzer"
    publishDir "${params.out_dir}/sam/cigar_analyzer", mode: 'copy'
    
    input:
        path sam


    output:
        path "*_top100.txt"                          , emit: comparison
        path "cigar_counts/combined_cigar_counts.txt", emit: stat 

    """
    python \\
        ${projectDir}/bin/cigar_analyzer.py \\
        --path ${params.testing_sam} \\
        --cigar ${params.cigar} \\
        --output . 
    head -15 cigar_counts/*.txt  >> cigar_counts/combined_cigar_counts.txt
    """
}