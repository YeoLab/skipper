import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact,chisquare
import math
import matplotlib.pyplot as plt
import seaborn as sns
def is_coding(s):
    if 'missense' in s or 'stop' in s or 'start' in s or 'frame' in s:
        return True
    else:
        return False

def plot_sorted(counts, ax, n=10, **kwargs):
    counts.loc[counts.sum(axis = 1).sort_values(ascending = False).index
    ].iloc[:n,:].plot.barh(stacked = True, **kwargs, ax = ax)

def plot_sorted_drop_column(counts, ax, n=10, columns_to_drop = ['Neutral'],**kwargs, ):
    subset = counts.loc[:,~counts.columns.isin(columns_to_drop)]
    subset.loc[subset.sum(axis = 1).sort_values(ascending = False).index,
    
    ].iloc[:n,:].plot.barh(stacked = True, **kwargs, ax = ax)

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

def test_enrichment_outliers(maf_vs_impact):
    ''' test enrichment to LoB and GoB variants'''
    stats = []
    background = 'neutral'
    for impact in ['GoB', 'LoB']:
        for index in maf_vs_impact.index:
            try:
                total = maf_vs_impact[background].sum()
                is_category = maf_vs_impact.loc[index, background]
        
                total_group = maf_vs_impact[impact].sum()
                is_category_group = maf_vs_impact.loc[index, impact]
                pv = pd.DataFrame(np.array([[is_category_group, total_group-is_category_group],
                               [is_category, total-is_category]]
                               ),
                                  index = [True, False],
                                  columns = [True, False],
                                 )
                if pv.ge(1).all().all():
                    pvalue, odds_ratio = testing(pv)
                    stats.append([index, impact, pvalue, odds_ratio])
            except Exception as e:
                print(e)
            
    stats = pd.DataFrame(stats, columns = ['MAF_bin', 'impact', 'pvalue', 'odds_ratio'])
    return stats
if __name__ == '__main__':
    clinvar_df = pd.read_csv(sys.argv[1])
    exp = Path(sys.argv[1]).name.split('.')[0]
    gnomad_output_dir = Path(sys.argv[2])
    gnomad = pd.read_csv(sys.argv[3], na_values = '.')
    clinvar_vep = pd.read_csv(sys.argv[4],
                              sep = '\t',
                      comment = '#',
                    names = ['Uploaded_variation','Location','Allele',
                      'Gene','Feature',
                      'Feature_type','Consequence','cDNA_position','CDS_position',
                      'Protein_position','Amino_acids','Codons','Existing_variation',
                      'Extra'])
    outdir = Path(sys.argv[5])

    # join gnomad
    alldf = []
    for f in list(gnomad_output_dir.glob('*.vcf')):
        
        chr_n = f.name.split('.')[0]
        print(chr_n)
        subset = clinvar_df[clinvar_df['CHROM']==f'chr{chr_n}']
        gdf = pd.read_csv(f, sep = '\t', na_values = '.',
                        names =['CHROM','POS','ID','REF','ALT','INFO/AC','INFO/AN','INFO/MR','INFO/AR','INFO/MG','INFO/MC']
            )
        alldf.append(subset.merge(gdf, left_on = ['CHROM','POS','REF','ALT'],
                    right_on = ['CHROM','POS','REF','ALT'],
                    how = 'left', suffixes = ('', '_gnomad')))
        
    clinvar_df = pd.concat(alldf, axis = 0)

    # get std from MAF>0.01
    std = gnomad.loc[(gnomad['INFO/AC']/gnomad['INFO/AN'])>0.01,'delta_score'].std()
    if math.isnan(std):
        std = gnomad['delta_score'].std()

    # join data
    clinvar_df['Consequence'] = clinvar_df['ID'].map(
        clinvar_vep.groupby(by = 'Uploaded_variation')['Consequence'].agg(
        lambda x: ','.join(list(set(','.join(x).split(',')))))
    )

    # clean up CLNSIG
    clinvar_df['INFO/CLNSIG']=clinvar_df['INFO/CLNSIG'].replace(
        'Likely_benign', 'Benign/Likely_benign').replace(
        'Benign', 'Benign/Likely_benign').replace(
        'Pathogenic', 'Pathogenic/Likely_pathogenic').replace(
        'Likely_pathogenic', 'Pathogenic/Likely_pathogenic').replace(
        'Conflicting_interpretations_of_pathogenicity', 'Uncertain_significance')
    clinvar_df['INFO/CLNSIG'] = pd.Categorical(clinvar_df['INFO/CLNSIG'], 
    categories = ['Benign/Likely_benign', 'Uncertain_significance', 'Pathogenic/Likely_pathogenic'],
    ordered = True)
    clinvar_df['CLNSIG_rank'] = clinvar_df['INFO/CLNSIG'].rank(method = 'dense')

    clinvar_df['INFO/CLNDN'] = clinvar_df['INFO/CLNDN'].str.replace('not_provided', 'not_specified')
    clinvar_df['impact'] = pd.cut(clinvar_df['delta_score'], 
                                    bins = [-np.inf,-2*std,2*std, np.inf],
                            labels = ['LoB', 'neutral', 'GoB']
                               )
    clinvar_df['delta_zscore'] = pd.cut(clinvar_df['delta_score'], 
                                    bins = [-np.inf,-2*std,-std, std ,2*std, np.inf],
                            labels = ['LoB', 'weak LoB', 'Neutral', 'weak GoB', 'GoB']
                               )
    clinvar_df['is_coding']=clinvar_df['Consequence'].fillna('').apply(lambda x: is_coding(x))
    clinvar_df['CLNDN_parsed'] = clinvar_df['INFO/CLNDN'].str.split('|').tolist()
    clinvar_df.to_csv(outdir / f'{exp}.clinvar_variants.csv')

    # plotting
    # perform counting
    clinvar_df_exploded = clinvar_df.explode('CLNDN_parsed')
    clinsig_count = pd.pivot_table(clinvar_df_exploded, index = 'CLNDN_parsed', 
            columns = ['is_coding','INFO/CLNSIG'],
           aggfunc = 'size')

    f, ax = plt.subplots(2,1, figsize = (6,6))
    try:
        plot_sorted(clinsig_count[True], ax = ax[0], cmap = 'coolwarm', legend = False)
    except KeyError as e:
        print(e)
    ax[0].set_title('coding')
    try:
        plot_sorted(clinsig_count[False], ax = ax[1], cmap = 'coolwarm')
    except KeyError as e:
        print(e)
    ax[1].set_title('non-coding')
    ax[1].set_xlabel('Number Variants')
    sns.despine()
    plt.savefig(outdir / f'{exp}.clinvar_CLINSIC_counts.pdf', bbox_inches='tight')
    clinsig_count.to_csv(outdir / f'{exp}.clinvar_CLINSIC_counts.csv')

    impact_count = pd.pivot_table(clinvar_df_exploded, index = 'CLNDN_parsed', 
            columns = ['is_coding','delta_zscore'],
           aggfunc = 'size')
    impact_count.to_csv(outdir / f'{exp}.clinvar_impact_counts.csv')
    f, ax = plt.subplots(2,1, figsize = (6,6))
    try:
        plot_sorted_drop_column(impact_count[True], ax = ax[0], cmap = 'coolwarm', legend = False)
    except KeyError as e:
        print(e)
    ax[0].set_title('coding')
    try:
        plot_sorted_drop_column(impact_count[False], ax = ax[1], cmap = 'coolwarm')
    except KeyError as e:
        print(e)
    ax[1].set_title('non-coding')
    ax[1].set_xlabel('Number Variants')
    sns.despine()
    plt.savefig(outdir / f'{exp}.clinvar_impact_counts.pdf', bbox_inches='tight')

    impact_vs_sig = pd.pivot_table(clinvar_df_exploded, columns = 'impact', 
            index = ['INFO/CLNSIG'],
           aggfunc = 'size')
    impact_enrichment = test_enrichment_outliers(impact_vs_sig)
    impact_enrichment.to_csv(outdir / f'{exp}.clinvar_impact_enrichment.csv')