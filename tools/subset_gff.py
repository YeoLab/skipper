# see https://bioinformatics.stackexchange.com/questions/13942/how-to-subset-genes-and-its-nested-features-from-a-gff-file-using-a-gene-list
import pandas as pd
import pyranges as pr
import argparse

parser = argparse.ArgumentParser(description='Filter GFF3 to relevant transcripts')
parser.add_argument('-a', '--full_annotation', metavar='gff3', required = True,
                    help='path to GFF3 annotation file')
parser.add_argument('-t', '--tpm_threshold', metavar="tpm", type=float, default = -1,
                    help='minimum TPM for GFF3 filtering')
parser.add_argument('-q', '--quant', metavar='quant.sf', required = True,
                    help='Salmon quantification file')
parser.add_argument('-o', '--subset_annotation', metavar='gff3', required = True,
                    help='path for subsetted GFF3 file')

args = parser.parse_args()

# full_annotation = "gencode.v38.annotation.gff3"
full_annotation = args.full_annotation
# tpm_threshold = 0.1
tpm_threshold = args.tpm_threshold
# quant = "quants/k562_totalrna/quant.sf"
quant = args.quant
# subset_annotation = "temp.gff3"
subset_annotation = args.subset_annotation

quant_data = pd.read_table(quant)
transcript_subset = set(quant_data[quant_data.TPM > tpm_threshold].Name)

gr_full = pr.read_gff3(full_annotation)

gr_subset = gr_full[pd.Series(transcript in transcript_subset for transcript in gr_full.transcript_id)]

gr_subset.to_gff3(subset_annotation)
