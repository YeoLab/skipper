# Load necessary packages 
import pandas as pd
import pyranges as pr
import argparse

# set up arguments. 
parser = argparse.ArgumentParser(description='Filter GFF3 to relevant transcripts')
parser.add_argument('-a', '--full_annotation', required = True,
                    help='path to GFF3 annotation file')
parser.add_argument('-t', '--tpm_threshold', type=float, default = 1,
                    help='minimum TPM for GFF3 filtering')
parser.add_argument('-s', '--source', required = True, choices=['salmon', 'GTEx', 'CCLE'],
                   help='source of TPM quantification, one of either "salmon", "GTEx", "CCLE"')
parser.add_argument('-c', '--cell_type', required=False,
                    help='Cell type from the GTEx or CCLE database (e.g. column) Required if source is GTEx or CCLE')
parser.add_argument('-q', '--quant', required = True,
                    help='TPM quantification file (.sf, .gtc, .txt)')
parser.add_argument('-o', '--subset_annotation', metavar='gff3', required = True,
                    help='path for subsetted GFF3 file')

# Load arguments. 
args = parser.parse_args()
full_annotation = args.full_annotation
tpm_threshold = args.tpm_threshold
source = args.source
cell_type = args.cell_type
quant = args.quant
subset_annotation = args.subset_annotation

# Conditional requirement. 
if args.source in {'GTEx', 'CCLE'} and not args.cell_type:
    parser.error("--cell_type is required when --source is GTEx or CCLE")

# Define helper functions. 
def load_expression_table(path):
    with open(path, 'r') as f:
        first_line = f.readline().strip()
    
    # Detect GCT
    if first_line.startswith("#"):
        df = pd.read_csv(path, sep='\t', skiprows=2)
    else:
        df = pd.read_csv(path, sep='\t')
    
    return df

def strip_after_period(values):
    return pd.Series(values).astype(str).str.split(".", n=1).str[0]


if source == "GTEx":
    # Load TPM data.
    quant_data = load_expression_table(quant)
    
    # Load the GFF file.
    gr_full = pr.read_gff3(full_annotation)
    
    # Strip version numbers from gene IDs in the expression table.
    quant_data["Name_stripped"] = strip_after_period(quant_data["Name"])
    
    # Find the subset of genes with a TPM above the threshold.
    gene_subset = set(quant_data.loc[quant_data[cell_type] > tpm_threshold, "Name_stripped"])
    
    # Strip version numbers from gene IDs in the GFF and subset.
    gff_gene_ids = strip_after_period(gr_full.gene_id)
    gr_subset = gr_full[gff_gene_ids.isin(gene_subset)]
    
    # Save the subset.
    gr_subset.to_gff3(subset_annotation)

if source == "CCLE":

    # Load TPM data.
    quant_data = load_expression_table(quant)
    
    # Load the GFF file.
    gr_full = pr.read_gff3(full_annotation)
    
    # Strip version numbers from gene IDs in the expression table.
    quant_data["Name_stripped"] = strip_after_period(quant_data["gene_id"])
    
    # Find the subset of genes with a TPM above the threshold.
    gene_subset = set(quant_data.loc[quant_data[cell_type] > tpm_threshold, "Name_stripped"])
    
    # Strip version numbers from gene IDs in the GFF and subset.
    gff_gene_ids = strip_after_period(gr_full.gene_id)
    gr_subset = gr_full[gff_gene_ids.isin(gene_subset)]
    
    # Save the subset.
    gr_subset.to_gff3(subset_annotation)

if source == "salmon":
    quant_data = pd.read_table(quant)
    transcript_subset = set(quant_data[quant_data.TPM > tpm_threshold].Name)
    
    gr_full = pr.read_gff3(full_annotation)
    
    gr_subset = gr_full[pd.Series(transcript in transcript_subset for transcript in gr_full.transcript_id)]
    
    gr_subset.to_gff3(subset_annotation)