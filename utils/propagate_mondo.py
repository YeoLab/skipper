from pathlib import Path
import numpy as np
from nxontology.imports import from_file
import seaborn as sns
import warnings
import networkx as nx
import pandas as pd
import sys
from scipy.stats import fisher_exact,chisquare
import numpy as np
warnings.filterwarnings("ignore")

def propagate_term(query):
    terms = set()
    upstream_edges = list(nx.edge_dfs(nxo.graph,query, orientation='reverse'))
    for e in upstream_edges:
        terms.add(e[0])
    
    return terms

def testing(df):
    ''' perform fisher exact or chisq given contingency table'''
    if df.shape != (2,2):
        # some outcomes are unobserved
        print(df.shape)
        print('no binding is observed in')
        return 1, np.nan
    if df.le(5).any().any():
        odds_ratio, pvalue = fisher_exact(df)
    else:
        chi, pvalue = chisquare(df.loc[True], (df.loc[True].sum())*df.loc[False].div(df.loc[False].sum()))
        odds_ratio = (df.loc[True, True]/df.loc[True, False])/(df.loc[False, True]/df.loc[False, False])
    return pvalue, odds_ratio

def propagate_counts(df, name):
    df['CLNDISDB_parsed'] = df['INFO/CLNDSDB'].str.replace('|', ',').str.split(',')
    df.reset_index(inplace = True)
    df.drop('index',axis = 1, inplace = True)
    
    df['variant_index']=df.index
    assert not df['variant_index'].duplicated().any()
    df['MAF']=(df['INFO/AC']/df['INFO/AN']).fillna(0)
    

    exploded_df = df.explode('CLNDISDB_parsed')
    exploded_df['CLNDISDB_parsed'] = exploded_df['CLNDISDB_parsed'].fillna('')
    subset = exploded_df.loc[exploded_df['CLNDISDB_parsed'].str.contains('MONDO')]
    need_mapping = exploded_df.loc[~exploded_df['CLNDISDB_parsed'].str.contains('MONDO')]
    mapped = need_mapping.merge(xref_mapping, left_on = 'CLNDISDB_parsed', right_on = 'xref')
    subset = pd.concat([subset.rename({'CLNDISDB_parsed':'MONDO'}, axis = 1), mapped], axis = 0, ignore_index = True)
    subset['MONDO_term'] = subset['MONDO'].str[6:]
    subset['terms']=subset['MONDO_term'].apply(propagate_term)
    variant_to_mondo_set = subset.groupby(by = ['variant_index'])['terms'].agg(lambda d: list(set().union(*d)))
    df['mondo_set'] = df['variant_index'].map(variant_to_mondo_set)
    
    try:
        cnts = pd.pivot_table(df.loc[(df['MAF']<0.01)&(~df['is_coding'])].explode('mondo_set'),
               index = 'mondo_set',
               columns = 'impact',
               aggfunc= 'size').fillna(0)
        
        
    except Exception as e:
        print(e)
    return df, cnts

def get_level(i):
    try:
        level = nx.shortest_path_length(nxo.graph, 'MONDO:0700096', i)
    except:
        level = -1
        
    return level
def enrichment_analysis(cnts):
    ''' for each level of the ontology, perform enrichment analysis '''
    stats = []
    for name, group in cnts.groupby(by = 'level'):
        total_in_level = group.sum()
        for index, row in group.iterrows():
            residual = total_in_level - row
            for impact in ['LoB', 'GoB']:
                if row[impact]>0:
                    pivot_table = pd.DataFrame(
                        [
                        [row[impact], row['neutral']],
                        [residual[impact], residual['neutral']]
                        ],
                        index = [True, False],
                        columns = [True, False]
                    )
                    pv, odds_ratio = testing(pivot_table)
                    stats.append([impact, index, pv, odds_ratio])
    stats = pd.DataFrame(stats, columns = ['impact', 'term', 'pvalue', 'odds_ratio'])
    stats = stats.merge(cnts, left_on = 'term', right_index = True)
    return stats

if __name__ == '__main__':
    nxo = from_file("https://purl.obolibrary.org/obo/mondo.owl")
    xref_mapping = pd.read_csv('/tscc/nfs/home/hsher/ontology/mondo/mondo_href.csv', index_col = 0)
    xref_mapping['xref'] = xref_mapping['xref'].str.replace('MESH:', 'MeSH:').str.replace('MEDGEN:', 'MedGen:')

    clinvar_variants = pd.read_csv(sys.argv[1], index_col = 0)
    name = Path(sys.argv[1]).name.split('.')[0]
    outprefix = sys.argv[2]

    df, cnts = propagate_counts(clinvar_variants, name)
    cnts['level']=pd.Series(cnts.index).apply(get_level).tolist()
    cnts=cnts.loc[cnts['level']>0]

    stats = enrichment_analysis(cnts)
    stats['exp'] = name

    print('saving to', f'{outprefix}.clinvar_mondo_tagged.csv')
    df.to_csv(f'{outprefix}.clinvar_mondo_tagged.csv')
    cnts.to_csv(f'{outprefix}.clinvar_mondo_cnts.csv')
    stats.to_csv(f'{outprefix}.clinvar_mondo_stats.csv')