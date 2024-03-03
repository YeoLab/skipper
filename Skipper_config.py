# Configuration file 

import os
import sys
import glob

# Adjust REPO_PATH if necessary to make sure it contains the path to the pulled skipper repository when starting a skipper run
REPO_PATH = "."

# Input data
########################################
# Files provided in the annotation folder, only adjust paths if needed 

# Use a GFF filtered for genes expressed in the cell type of interest
GFF = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.gff3.gz"
PARTITION = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.bed.gz"
# Skipper will partition the transcriptome and create feature annotations from the GFF
FEATURE_ANNOTATIONS = REPO_PATH + "/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.features.tsv.gz"
# Downloaded from UCSC genome browser
CHROM_SIZES = REPO_PATH + "/annotations/hg38.chrom.sizes"
# Customizable, with defaults
BLACKLIST = REPO_PATH + "/annotations/encode3_eclip_blacklist.bed" # set to None for no blacklisting
GENE_SETS = REPO_PATH + "/annotations/c5.go.v7.5.1.symbols.gmt"
GENE_SET_REFERENCE = REPO_PATH + "/annotations/encode3_go_terms.reference.tsv.gz"
GENE_SET_DISTANCE = REPO_PATH + "/annotations/encode3_go_terms.jaccard_index.rds"
REPEAT_BED = REPO_PATH + "/annotations/repeatmasker.grch38.sort.unique.bed.gz"
# Ranked list of gene and transcript types found in GFF annotations
ACCESSION_RANKINGS = REPO_PATH + "/annotations/accession_type_ranking.txt"
TOOL_DIR = REPO_PATH + "/tools"

########################################
# Files not provided, setup required 

MANIFEST = REPO_PATH + "/example/Example_manifest.csv"
# General resources
REPEAT_TABLE = REPO_PATH + "/annotations/repeatmasker.grch38.tsv.gz"
# indexed genome fasta
GENOME = REPO_PATH + "/annotations/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
# STAR reference
STAR_DIR = REPO_PATH + "/annotations/genome_ref"
# Use singularity, if False, make sure to download UMICollapse v1.0.0 into installation folder
SINGULARITY = False

########################################
# Customizable Parameters

# Information about CLIP library
UMI_SIZE = 10
# Single-end: enter 1. Paired-end: enter read (1 or 2) corresponding to crosslink site
INFORMATIVE_READ = 1
# Internal use
UNINFORMATIVE_READ = 3 - INFORMATIVE_READ
# Use multiple input replicates to estimate overdispersion (preferred), or use multiple CLIP replicates
# Skipper requires replicates to model the variance in read counts.
OVERDISPERSION_MODE = "input" # input or clip