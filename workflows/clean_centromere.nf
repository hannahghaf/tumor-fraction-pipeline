nextflow.enable.dsl=2

workflow clean_centromere {
    take:
    raw_centromere_ch

    main:
    centromere_bed_ch = clean_cent(raw_centromere_ch)

    emit:
    centromere_bed_ch
}

process clean_cent {
    input:
    path centromere_file

    output:
    path "centromere.hg38.bed"

    script:
    """
    # remove header
    tail -n +2 ${centromere_file} | cut -f1-3 > centromere.hg38.bed
    """
}