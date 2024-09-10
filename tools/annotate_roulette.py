import sys
import pandas as pd
from pathlib import Path
import numpy as np
def sigmoid(z):
    return 1/(1 + np.exp(-z))

def join_data(indir, exp, chrom):
    ref = pd.read_csv(indir / 'output/variants' / 'roulette' / f'{exp}.{chrom}.ref.score.txt', sep = '\t', names = ['ID', 'score'])
    ref.rename({'score': 'score_ref'}, axis = 1, inplace = True)
    ref.drop('ID', axis = 1, inplace = True)
    alt = pd.read_csv(indir / 'output/variants' / 'roulette' / f'{exp}.{chrom}.alt.score.txt', sep = '\t', names = ['ID', 'score'])
    alt.rename({'score': 'score_alt'}, axis = 1, inplace = True)
    alt.drop('ID', axis = 1, inplace = True)
    annot = pd.read_csv(indir / 'output/variants' / 'roulette' / f'{exp}.{chrom}.csv',
                        index_col = 0
                       ).rename({'5': 'FILTER',
                                 '6': 'INFO/MR',
                                 '7': 'INFO/AR',
                                 '8': 'INFO/MG',
                                 '9': 'INFO/MC'
                                }, axis = 1
                               )
    print(chrom)
    assert ref.shape[0]==alt.shape[0]
    assert annot.shape[0]==ref.shape[0]
    
    # calculate MAF and binning
    df = pd.concat([ref,alt, annot], axis = 1)
    
    df['score_alt'] = sigmoid(df['score_alt'])
    df['score_ref'] = sigmoid(df['score_ref'])
    df['delta_score']=df['score_alt']-df['score_ref']

    return df

if __name__ =='__main__':

    indir = Path(sys.argv[1])
    exp = sys.argv[2]

    # read all roulette data 
    alldf = []
    for n in range(1,23):
        try:
            alldf.append(join_data(indir, exp, f'chr{n}'))
        except Exception as e:
            print(e, 'chr{n}')
    df = pd.concat(alldf, axis = 0)

    # merge with annotation
    site_annotation = pd.read_csv(indir / 'output/finemapping/mapped_sites/' / f'{exp}.finemapped_windows.annotated.tsv',
                             sep = '\t')
    
    df = df.merge(site_annotation[['name', 'feature_type_top',
       'feature_types', 'gene_name', 'gene_id', 'transcript_ids',
       'gene_type_top', 'transcript_type_top', 'gene_types',
       'transcript_types']], left_on = 'name', right_on = 'name')
    
    # add gnomAD annotation
    gnomad = pd.read_csv(indir / 'output/variants' / 'gnomAD' / f'{exp}.csv',
                    index_col = 0).rename({'5': 'INFO/AC', '6':'INFO/AN'}, axis = 1)
    gnomad['MAF']=gnomad['INFO/AC']/gnomad['INFO/AN']
    gnomad['MAF_bin'] = pd.cut(gnomad['MAF'], bins = [0, 1e-3, 1e-2, 1])
    bins = gnomad['MAF_bin'].cat.categories.tolist()
    bins = ['singleton']+bins
    gnomad['MAF_bin']=pd.Categorical(gnomad['MAF_bin'], categories = bins, ordered = True)
    gnomad.loc[gnomad['INFO/AC']==1, 'MAF_bin']='singleton'
    gnomad['MAF_bin_rank']=gnomad['MAF_bin'].rank(method = 'dense')
    df = df.merge(gnomad[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'MAF', 'MAF_bin', 'MAF_bin_rank',
                          'INFO/AC', 'INFO/AN']], how = 'left', 
         left_on = ['CHROM', 'POS', 'REF', 'ALT'],
         right_on = ['CHROM', 'POS', 'REF', 'ALT']
        )
    df['MAF_bin']=pd.Categorical(df['MAF_bin'], categories = ['unobserved']+bins, ordered = True)
    df['MAF_bin'].fillna('unobserved', inplace = True)

    # add gene-wise election constrain
    gnomad_constrain = pd.read_csv('/tscc/projects/ps-yeolab5/hsher/gnomAD/gnomad.v4.0.constraint_metrics.tsv', sep = '\t')
    gnomad_constrain = gnomad_constrain.loc[(gnomad_constrain['mane_select'])&(gnomad_constrain['transcript'].str.contains('ENST'))]
    gnomad_constrain['lof.pLI bins']=pd.cut(gnomad_constrain['lof.pLI'], bins = [0,0.1,0.9,1],
                            labels = ['LoF tolerant', 'Ambiguous', 'LoF intolerant'])
    df['lof.pLI bins']=df['gene_name'].map(gnomad_constrain.dropna().set_index('gene')['lof.pLI bins'])
    df['mis_pphen.oe']=df['gene_name'].map(gnomad_constrain.dropna().set_index('gene')['mis_pphen.oe'])
    df['syn.oe']=df['gene_name'].map(gnomad_constrain.dropna().set_index('gene')['syn.oe'])

    # filter for good quality stuffs
    filtered_df = df.loc[df['FILTER']!='low']
    
    df.to_csv(indir / 'output/variants/roulette' / f'{exp}.full.csv.gz',compression='gzip')
    filtered_df.to_csv(indir / 'output/variants/roulette' / f'{exp}.filtered.csv.gz',compression='gzip')





