nextflow.enable.dsl=2

workflow run_ichorCNA {
    take:
    ichor_input_ch
    
    main:
    ichor_output_ch = run_ichor(ichor_input_ch)
    tf_files_ch = extract_tumor_fraction(ichor_output_ch)
    tf_table_ch = merge_tumor_fractions(tf_files_ch.collect())

    emit:
    tumor_fraction_table = tf_table_ch
}

process run_ichor {
    tag "$wig_file.simpleName"

    input:
    tuple path(wig_file), path(centromere_bed), path(gc_wig), path(map_wig), path(pon_rds)

    output:
    path "${wig_file.simpleName}_ichor_results"

    script:
    """
    mkdir -p ${wig_file.simpleName}_ichor_results

    Rscript /opt/ichorCNA/scripts/runIchorCNA.R \
        --id ${wig_file.simpleName} \
        --WIG ${wig_file} \
        --gcWig ${gc_wig} \
        --mapWig ${map_wig} \
        --centromere ${centromere_bed} \
        --normalPanel ${pon_rds} \
        --outDir ${wig_file.simpleName}_ichor_results \
        --includeHOMD FALSE \
        --ploidy 2 \
        --maxCN 5 \
        --chrTrain "c(1:22)" \
        --estimateNormal TRUE \
        --estimatePloidy TRUE \
        --estimateScPrevalence TRUE \
        --txnE 0.9999 \
        --txnStrength 10000 \
        --normal 0.5 \
        --minSegmentBins 50 \
        --maxFracCNASubclone 1 \
        --plotFileType png
    """

}

process extract_tumor_fraction {
    tag "$ichor_out_dir.simpleName"

    input:
    path ichor_out_dir

    output:
    path "*.tumor_fraction_summary.tsv"

    script:
    """
    sample_name=\$(basename ${ichor_out_dir} | sed 's/_ichor_results//')
    tf=\$(grep 'Tumor Fraction:' ${ichor_out_dir}/*.params.txt | awk '{print \$NF}')
    echo -e "sample\\ttumor_fraction" > \${sample_name}.tumor_fraction_summary.tsv
    echo -e "\$sample_name\\t\$tf" >> \${sample_name}.tumor_fraction_summary.tsv
    """
}

process merge_tumor_fractions {
    tag "merge tumor fractions"
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path tf_files

    output:
    path "tumor_fraction_summary.tsv"

    script:
    """
    head -n 1 \$(ls ${tf_files} | head -n 1) > tumor_fraction_summary.tsv
    tail -n +2 -q ${tf_files} >> tumor_fraction_summary.tsv
    """
}