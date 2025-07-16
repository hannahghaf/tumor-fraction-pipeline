# Tumor Fraction Estimation Pipeline

This repository contains a Nextflow pipeline that estimates tumor fraction from BAM files using the [ichorCNA] (https://github.com/broadinstitute/ichorCNA) algorithm.
The pipeline processes raw BAM (HG19) files, converts them to BED (both HG19 and HG38) and WIG formats, and runs ichorCNA to output tumor fraction summaries.

## Setup

### Requirements
- [Nextflow] (https://www.nextflow.io/docs/latest/install.html)
- [Docker] (https://docs.docker.com/get-docker/)
- Reference files (see below)

### Input files
- `*.bam` files (in `params.input_dir`)
- `chain_file` for LiftOver (e.g., `hg19ToHg38.over.chain.gz`)
- `chrom_sizes` file (e.g., `hg38.chrom.sizes`)
- `centromere`, `gcWig`, `mapWig`, and `normalPanel.rds` files for ichorCNA

## Usage

### Basic command
```bash
nextflow run main.nf -profile standard \
	--input_dir data/bams \
	--chain_file data/hg19ToHg38.over.chain.gz \
	--chrom_sizes data/hg38.chrom.sizes \
	--raw_centromere data/centromere_hg38.txt \
	--gc_wig data/gc.wig \
	--map_wig data/map.wig \
	--pon_rds data/normalPanel.rds \
	--out_dir results
```

## Output
- `results/tumor_fraction_summary.tsv` -> table of tumor fractions
- `bed/` `bed_hg38/` `wig/` files in `results/`
- intermediate files in `.nextflow/` and `work/`
