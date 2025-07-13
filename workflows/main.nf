nextflow.enable.dsl=2

params.input_dir = "./data/bam/200k"
params.out_dir = "./results"

workflow {
    bam_files = Channel.fromPath("${params.input_dir}/*.bam")

    convert_bam_to_bed(bam_files)
}

process convert_bam_to_bed {
    container 'my-bioimg:latest' //docker container w bedtools to convert bam to bed

    errorStrategy = 'retry'
    maxRetries = 2

    tag "$bam"
    publishDir "${params.out_dir}/bed", mode: 'copy'

    input:
    path bam

    output:
    path "*.bed"
    path "exit_status.txt"

    script:
    """
    bedtools bamtobed -i ${bam} > ${bam.simpleName}.bed

    echo \$? > exit_status.txt
    """
}