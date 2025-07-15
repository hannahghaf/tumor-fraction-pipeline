nextflow.enable.dsl=2

include { bamtobed } from './workflows/bamtobed.nf'
include { bedtowig } from './workflows/bedtowig.nf'

workflow {
    // get input BAM files from params.input_dir
    def bam_files_ch = Channel.fromPath("${params.input_dir}/*.bam")
    def chain_file_ch = Channel.value(file(params.chain_file))

    // bam to bed workflow BAM -> BED (HG38)
   def hg38_bed_files_ch = bamtobed(bam_files_ch, chain_file_ch)

   // bed to wig workflow for ichorCNA input prep
   def chrom_sizes_ch = Channel.fromPath(params.chrom_sizes)
   def wig_file_ch = bedtowig(hg38_bed_files_ch, chrom_sizes_ch)
}

workflow.onComplete {
    if (workflow.success) {
        println "Pipeline completed successfully at ${workflow.complete}"
    }
}