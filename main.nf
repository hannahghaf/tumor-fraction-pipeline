nextflow.enable.dsl=2

include { bamtobed } from './workflows/bamtobed.nf'
include { bedtowig } from './workflows/bedtowig.nf'
include { run_ichorCNA } from './workflows/run_ichorCNA.nf'
include { clean_centromere } from './workflows/clean_centromere.nf'

workflow {
    // get input BAM files from params.input_dir
    def bam_files_ch = Channel.fromPath("${params.input_dir}/*.bam")
    def chain_file_ch = Channel.value(file(params.chain_file))

    // bam to bed workflow BAM -> BED (HG38)
   def hg38_bed_files_ch = bamtobed(bam_files_ch, chain_file_ch)

   // bed to wig workflow for ichorCNA input prep
   def chrom_sizes_ch = Channel.value(file(params.chrom_sizes))
   def wig_file_ch = bedtowig(hg38_bed_files_ch, chrom_sizes_ch)

   // clean centromere file for ichor input
   def raw_centromere_ch = Channel.fromPath(params.raw_centromere)
   def centromere_bed_ch = clean_centromere(raw_centromere_ch)

   def gc_wig_ch = Channel.fromPath(params.gc_wig)
   def map_wig_ch = Channel.fromPath(params.map_wig)
   def pon_rds_ch = Channel.fromPath(params.pon_rds)

   // run ichorCNA and output results
   def ichor_input_ch = wig_file_ch
    .combine(centromere_bed_ch)
    .map { wig, centromere ->
        tuple(wig, centromere, file(params.gc_wig), file(params.map_wig), file(params.pon_rds))
    }
   def ichorCNA_output = run_ichorCNA(ichor_input_ch)
}

workflow.onComplete {
    if (workflow.success) {
        println "Pipeline completed successfully at ${workflow.complete}"
    }
}