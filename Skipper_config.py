# Configuration file 

import os
import sys
import glob

REPO_PATH = "/tscc/projects/ps-yeolab4/software/skipper/d0055ff/bin/skipper"
# Input data
# Use a GFF filtered for genes expressed in the cell type of interest
# Skipper will partition the transcriptome and create feature annotations from the GFF
MANIFEST = "manifest.csv"
GFF = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.gff3.gz"
PARTITION = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.bed.gz"
FEATURE_ANNOTATIONS = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.features.tsv.gz"

# Information about CLIP library
UMI_SIZE = 10
# Single-end: enter 1. Paired-end: enter read (1 or 2) corresponding to crosslink site
INFORMATIVE_READ = 1
# Use multiple input replicates to estimate overdispersion (preferred), or use multiple CLIP replicates
# Skipper requires replicates to model the variance in read counts.
OVERDISPERSION_MODE = "input" # input or clip

# Paths to custom tools
# Path to tools packaged with Skipper
TOOL_DIR = REPO_PATH + "/tools"
# Path to additional installed tools (enter "." for the current directory if all tools are globally installed)
EXE_DIR = ""

# General resources
# STAR reference
STAR_DIR = "/tscc/projects/ps-yeolab3/bay001/annotations/GRCh38/star_2_7_gencode40_sjdb"
# Path to R with installed dependencies: VGAM, viridis, fgsea, GenomicRanges, ggrepel, RColorBrewer, Rtsne, ggupset, cowplot, ggdendro
R_EXE = ""
# Downloaded from https://github.com/Daniel-Liu-c0deb0t/UMICollapse
UMICOLLAPSE_DIR = ""
# Java executable for UMICollapse: enter "java" for globally installed version
JAVA_EXE = ""
# indexed genome fasta
GENOME = "/tscc/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
# Downloaded from UCSC genome browser
CHROM_SIZES = "/tscc/projects/ps-yeolab3/bay001/annotations/GRCh38/star_2_7_gencode40_sjdb/chrNameLength.txt"
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1435901349_MASaYmzfsYNwQY34JLgNir7HOyVh&clade=mammal&org=Human&db=hg38
REPEAT_TABLE = "/tscc/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/repeatmasker.grch38.tsv.gz"

# Customizable, with defaults
BLACKLIST = REPO_PATH + "/annotations/encode3_eclip_blacklist.bed" # set to None for no blacklisting
GENE_SETS = REPO_PATH + "/annotations/c5.go.v7.5.1.symbols.gmt"
GENE_SET_REFERENCE = REPO_PATH + "/annotations/encode3_go_terms.reference.tsv.gz"
GENE_SET_DISTANCE = REPO_PATH + "/annotations/encode3_go_terms.jaccard_index.rds"
REPEAT_BED = "/tscc/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/repeatmasker.grch38.sort.unique.bed.gz"
# Ranked list of gene and transcript types found in GFF annotations
ACCESSION_RANKINGS = REPO_PATH + "/annotations/accession_type_ranking.txt"

# Internal use
UNINFORMATIVE_READ = 3 - INFORMATIVE_READ
