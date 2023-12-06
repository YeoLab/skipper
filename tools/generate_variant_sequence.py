import pandas as pd
from pybedtools import BedTool
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def generate_variant_sequence(row):
    ref = row['REF']
    seq = row['seq']
    
    
    
    start = row['POS']-row['start']-1
    end = row['POS']-row['start']-1+len(ref)
    
    to_replace = seq[start:end]
    if to_replace != ref:
        print(to_replace, 'ref=', ref, start, end, len(seq))
        return None
    new_seq = seq[:start]+row['ALT']+seq[end:]
    return new_seq

def reverse_complement(string):
    newstr = ''
    mapper={'A':'T',
            'T':'A',
            'C': 'G',
            'G':'C',
           '.': '',
           'N': 'N'}
    for s in string[::-1]:
        newstr+=mapper[s]
    return newstr

def generate_variant_sequence_neg(row):
    ref = row['REF']
    seq = row['seq']
    
    
    
    start = row['end']-row['POS']
    end = row['end']-row['POS']-len(ref)
    
    to_replace = reverse_complement(seq[end+1:start+1])
    if to_replace != ref:
        print(to_replace, 'ref=',ref, start, end, len(seq))
        return None
    new_seq = seq[:end+1]+reverse_complement(row['ALT'])+seq[start+1:]
    return new_seq



def make_SeqRecord(df, seq_col):
    ''' convert df to a list of SeqRecord '''
    records = []
    for index, row in df.iterrows():
        records.append(
            SeqRecord(
                Seq(row[seq_col]),
                id=str(row['ID'])
                )
        )
    return records

if __name__ == '__main__':
    vcf_file = sys.argv[1]
    fa = sys.argv[2]
    bed = sys.argv[3]
    out_prefix = sys.argv[4]

    # read vcf
    variants_df = pd.read_csv(vcf_file,
        sep = '\t',
        header = None).drop_duplicates(subset = [0,1,2,3,4])
    variants = variants_df.iloc[:, :5]
    variants.columns = ['CHROM','POS','ID','REF', 'ALT']

    sequences = {}
    for record in SeqIO.parse(fa, format = 'fasta'):
        sequences[record.id]=str(record.seq)

    # get sequence of binding sites
    window_bed = BedTool(bed)
    window_df = window_bed.to_dataframe()
    window_df['seq']=sequences.values()
    window_bed = BedTool.from_dataframe(window_df)


if variants.empty:
    variant_df = pd.DataFrame(columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'name', 'variant_seq','feature_type_top', 'feature_types', 'gene_name',
    'transcript_types', 'transcript_type_top', 'INFO/AC', 'INFO/AN'])
    variant_df.to_csv(outf)
else:
    # finding window name
    variants['POS-1']=variants['POS']-1
    
    # get variants in window
    df = BedTool.from_dataframe(variants[['CHROM','POS-1','POS','ID','REF', 'ALT']]).intersect(
        window_bed, wb = True).to_dataframe(names = ['CHROM','POS-1', 'POS','ID','REF', 'ALT']+window_df.columns.tolist())
    
    df.dropna(subset = ['seq'],inplace = True)
    df=df.loc[df['REF'].str.len()<10]

    pos = df.loc[df['strand']=='+']
    neg = df.loc[df['strand']=='-']

    # generate variant sequence
    if not pos.empty:
        pos['variant_seq']=pos.apply(generate_variant_sequence, axis = 1)
    else:
        pos['variant_seq']=None
    if not neg.empty:
        neg['variant_seq']=neg.apply(generate_variant_sequence_neg, axis = 1)
    else:
        neg['variant_seq']=None

    cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'seq', 'variant_seq','name']
    variant_seq_df = pd.concat([pos[cols],
        neg[cols]],
        axis = 0)
    variant_seq_df.dropna(subset = ['variant_seq'], inplace = True)
    variant_seq_df = variant_seq_df.merge(variants_df, 
                                          left_on = ['CHROM', 'POS', 'ID', 'REF', 'ALT'],
                                          right_on = [0,1,2,3,4]).drop([0,1,2,3,4], axis = 1)
    
    variant_seq_df.to_csv(f'{out_prefix}.csv')

    # make fasta
    ref_seq = make_SeqRecord(variant_seq_df, 'seq')
    alt_seq = make_SeqRecord(variant_seq_df, 'variant_seq')

    SeqIO.write(ref_seq, f'{out_prefix}.ref.fa', format = 'fasta')
    SeqIO.write(alt_seq, f'{out_prefix}.alt.fa', format = 'fasta')
