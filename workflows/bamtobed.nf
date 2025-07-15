nextflow.enable.dsl=2

workflow bamtobed {
    take:
    bam_ch
    chain_file_ch

    main:
    // BAM -> BED
    bed_ch = convert_bam_to_bed(bam_ch)
    // BED (HG19) -> BED (HG38)
    hg38_bed_files_ch = liftover_bed_to_hg38(bed_ch, chain_file_ch)

    emit:
    hg38_bed_files_ch
}

process convert_bam_to_bed {
    container 'my-bioimg:latest' //docker container w bedtools to convert bam to bed

    errorStrategy = 'retry'
    maxRetries = 2

    tag "$bam" //log readability
    publishDir "${params.out_dir}/bed", mode: 'copy'

    input:
    path bam

    output:
    path "*.bed"

    script:
    """
    # HG19 BAM to HG19 BED
    bedtools bamtobed -i ${bam} > ${bam.simpleName}.bed 
    """
}

process liftover_bed_to_hg38 {
    container 'my-bioimg:latest'
    errorStrategy = 'retry'
    maxRetries = 2

    tag "$bed_file"
    publishDir "${params.out_dir}/bed_hg38", mode: 'copy'

    input:
    path bed_file
    path chain_file

    output:
    path "${bed_file.simpleName}.hg38.bed"

    script:
    """
    # liftOver with chain_file to map HG19 to HG38 coords
    liftOver ${bed_file} ${chain_file} ${bed_file.simpleName}.hg38.bed unmapped.bed
    """
}