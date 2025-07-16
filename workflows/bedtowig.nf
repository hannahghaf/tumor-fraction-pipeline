nextflow.enable.dsl=2

workflow bedtowig {
    take:
    bed_files_ch    // HG38 BED files
    chrom_sizes_ch  // chrom.sizes reference

    main:
    bins_bed_ch = make_bins_bed(chrom_sizes_ch)
    paired_inputs_ch = bed_files_ch.combine(bins_bed_ch)
    counts_bed_ch = count_reads_per_bin(paired_inputs_ch)
    wig_ch = convert_to_wig(counts_bed_ch)
    clean_wigs_ch = filter_wig_for_ichor(wig_ch)

    emit:
    wig_file = clean_wigs_ch
}

process make_bins_bed {
    // generate bins BED file

    tag "$chrom_sizes.simpleName"

    input:
    path chrom_sizes

    output:
    path "bins_1mb.bed"

    script:
    """
    # make BED file of 1Mb bins
    bedtools makewindows -g ${chrom_sizes} -w 1000000 > bins_1mb.bed
    """
}

process count_reads_per_bin {
    // count reads per bin

    tag "$bed_file.simpleName"

    input:
    tuple path(bed_file), path(bins_bed)

    output:
    path "${bed_file.simpleName}.counts.bed"

    script:
    """
    # count how many of BED intervals falls into each 1Mb bin
    bedtools intersect -a ${bins_bed} -b ${bed_file} -c > ${bed_file.simpleName}.counts.bed
    """
}

process convert_to_wig {
    // convert count BED to WIG

    tag "$counts_bed.simpleName"
    publishDir "${params.out_dir}/wig", mode: 'copy'

    input:
    path counts_bed

    output:
    path "${counts_bed.simpleName}.wig"

    script:
    """
    # turn counts into .wig format
    awk 'BEGIN {OFS="\\t"} {print "fixedStep chrom=" \$1 " start=" \$2+1 " step=1000000 span=1000000\\n" \$4}' ${counts_bed} > ${counts_bed.simpleName}.wig
    """
}

process filter_wig_for_ichor {
    tag "$wig_file.simpleName"
    
    input:
    path wig_file
    
    output:
    path "${wig_file.simpleName}.clean.wig"
    
    script:
    """
    # remove problematic chromosomes from WIG file
    awk '
    /^fixedStep/ {
        # filter to only autosomes and chrX/Y, exclude random/alt/unplaced
        if (\$0 ~ /chrom=chr([1-9]|1[0-9]|2[0-2]|X|Y)[ \t]/ && \$0 !~ /_random|_alt|_fix|_hap|chrUn_|chrEBV/) {
            print \$0
            keep = 1
        } else {
            keep = 0
        }
        next
    }
    # only print data lines for chrs to keep
    keep == 1 { print \$0 }
    ' ${wig_file} > ${wig_file.simpleName}.clean.wig
    
    echo "Original WIG lines: \$(wc -l < ${wig_file})"
    echo "Filtered WIG lines: \$(wc -l < ${wig_file.simpleName}.clean.wig)"
    """
}