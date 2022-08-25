# Configuration file 

import os
import sys
import glob

# Input data
# Use a GFF filtered for genes expressed in the cell type of interest
# Skipper will partition the transcriptome and create feature annotations from the GFF
MANIFEST = '/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/20220803_ribo/ribo_manifest_20220124.csv'
GFF = '/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/gencode.v41.annotation.gff3.gz'
PARTITION = '/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/gencode.v41.annotation.tiled_partition.bed.gz'
FEATURE_ANNOTATIONS = '/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/gencode.v41.annotation.tiled_partition.features.tsv.gz'
ACCESSION_RANKINGS = '/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/accession_type_ranking.txt'

# Information about CLIP library
UMI_SIZE = 10
# Single-end: enter 1. Paired-end: enter read (1 or 2) corresponding to crosslink site
INFORMATIVE_READ = 1
# Use multiple input replicates to estimate overdispersion (preferred), or use multiple CLIP replicates
# Skipper requires replicates to model the variance in read counts.
OVERDISPERSION_MODE = "input" # input or clip

# Paths to custom tools
CONDA_DIR = '~/miniconda3'
TOOL_DIR = '/projects/ps-yeolab3/eboyle/encode/pipeline/06_20220816/tools'
EXE_DIR = '/projects/ps-yeolab3/eboyle/software/bin/'

# General resources
# STAR reference
STAR_DIR = '/projects/ps-yeolab4/software/eclip/0.7.0/examples/inputs/star_2_7_gencode29_sjdb'
# Path to R with installed dependencies: VGAM, viridis, fgsea, GenomicRanges, ggrepel, RColorBrewer, Rtsne, ggupset
R_EXE = '/projects/ps-yeolab4/software/R-4.1.2/bin/Rscript'
# Downloaded from https://github.com/Daniel-Liu-c0deb0t/UMICollapse
UMICOLLAPSE_DIR = "/projects/ps-yeolab3/eboyle/software/UMICollapse"
# Java executable for UMICollapse
JAVA_EXE = '/projects/ps-yeolab3/eboyle/software/jdk-18.0.1.1/bin/java'
# indexed genome fasta
GENOME = '/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta'
# Downloaded from UCSC genome browser
CHROM_SIZES = '/projects/ps-yeolab3/eboyle/resources/hg38.chrom.sizes'
REPEAT_TABLE = '/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/repeatmasker.grch38.tsv.gz'
# Customizable, with defaults
BLACKLIST = '/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/20220728_encode3/encode3_eclip_blacklist.bed'
GENE_SETS = "/projects/ps-yeolab3/eboyle/encode/pipeline/06_20220816/tools/c5.go.v7.5.1.symbols.gmt"
GENE_SET_REFERENCE = '/projects/ps-yeolab3/eboyle/encode/pipeline/06_20220816/tools/encode3_go_terms.reference.tsv.gz'
GENE_SET_DISTANCE = '/projects/ps-yeolab3/eboyle/encode/pipeline/06_20220816/tools/encode3_go_terms.jaccard_index.rds'

# Internal use
UNINFORMATIVE_READ = 3 - INFORMATIVE_READ