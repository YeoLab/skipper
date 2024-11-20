import pandas as pd
from pathlib import Path
from pybedtools import BedTool
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import fisher_exact
from scipy.stats import chisquare

def sigmoid(z):
    return 1/(1 + np.exp(-z))

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

if __name__ == '__main__':
    indir = Path(sys.argv[1])
    exp = sys.argv[2]
    finemap_annotation = indir / 'output/finemapping/mapped_sites/' / f'{exp}.finemapped_windows.annotated.tsv'

    # load all files
    ref = pd.read_csv(indir / 'output/variants' / 'clinvar' / f'{exp}.ref.score.txt', sep = '\t', names = ['ID', 'score'])
    ref.rename({'score': 'score_ref'}, axis = 1, inplace = True)
    ref.drop('ID', axis = 1, inplace = True)
    alt = pd.read_csv(indir / 'output/variants' / 'clinvar' / f'{exp}.alt.score.txt', sep = '\t', names = ['ID', 'score'])
    alt.rename({'score': 'score_alt'}, axis = 1, inplace = True)
    alt.drop('ID', axis = 1, inplace = True)
    annot = pd.read_csv(indir / 'output/variants' / 'clinvar' / f'{exp}.csv', index_col = 0)
    site_annotation_uniq = pd.read_csv(finemap_annotation, sep = '\t')

    # load gnomAD test
    gnomAD_test = pd.read_csv(indir / 'output/variants/gnomAD_analysis' / f'{exp}.test_stats.csv', index_col = 0)
    # calculate  binning
    df = pd.concat([ref,alt, annot], axis = 1)
    df.rename({'5': 'INFO/CLNDN', '6': 'INFO/CLNVC', '7': 'INFO/CLNSIG',
          '8': 'INFO/CLNDISDB', '9': 'INFO/AF_ESP', '10': 'INFO/AF_EXAC',
          '11': 'INFO/AF_TGP', '12': 'INFO/ALLELEID'}
          , axis = 1, inplace = True)
    
    # clean up CLNSIG
    df['INFO/CLNSIG']=df['INFO/CLNSIG'].replace(
        'Likely_benign', 'Benign/Likely_benign').replace(
        'Benign', 'Benign/Likely_benign').replace(
        'Pathogenic', 'Pathogenic/Likely_pathogenic').replace(
        'Likely_pathogenic', 'Pathogenic/Likely_pathogenic').replace(
        'Conflicting_interpretations_of_pathogenicity', 'Uncertain_significance')
    df['INFO/CLNSIG'] = pd.Categorical(df['INFO/CLNSIG'], 
    categories = ['Benign/Likely_benign', 'Uncertain_significance', 'Pathogenic/Likely_pathogenic'],
    ordered = True)
    df['CLNSIG_rank'] = df['INFO/CLNSIG'].rank(method = 'dense')
    df['score_alt'] = sigmoid(df['score_alt'])
    df['score_ref'] = sigmoid(df['score_ref'])
    df['delta_score']=df['score_alt']-df['score_ref']

    
    # now annotate, all variants should be derived from finemapped windows
    assert df['name'].isin(site_annotation_uniq['name']).all()
    df['feature_type_top'] = df['name'].map(site_annotation_uniq.set_index('name')['feature_type_top'])
    df['transcript_type_top'] = df['name'].map(site_annotation_uniq.set_index('name')['transcript_type_top'])
    df['gene_name']=df['name'].map(site_annotation_uniq.set_index('name')['gene_name'])
    df['strand'] = df['name'].map(site_annotation_uniq.set_index('name')['strand'])

    # annotate variant type
    df.loc[(df['REF'].str.len()==1)&(df['ALT'].str.len()==1), 'TYPE']='SNV'
    df.loc[(df['REF'].str.len()!=1)|(df['ALT'].str.len()!=1), 'TYPE']='INDEL'

    # annotat impact
    # annotate LoF and GoF based on gnomAD std
    with open(indir / 'output/variants' / 'gnomAD' / f'{exp}.variants_scores.tsv', 'r') as f:
        for line in f:
            std=float(line.split(':')[1].rstrip())
            break
    df.loc[df['delta_score']<-2*std, 'impact'] = 'LoF'
    df.loc[df['delta_score']>2*std, 'impact'] = 'GoF'
    df['impact'].fillna('neutral', inplace = True)

    # Summary statistics: Distribution of delta
    with open(indir / 'output/variants' / 'gnomAD' / f'{exp}.variants_scores.tsv', 'r') as f:
        for line in f:
            if 'std' in line:
                std=float(line.split(':')[1].rstrip())
            elif 'mean' in line:
                mean=float(line.split(':')[1].rstrip())
            else:
                break
    
    df.loc[df['delta_score']<-2*std, 'impact'] = 'LoF'
    df.loc[df['delta_score']>2*std, 'impact'] = 'GoF'
    df['impact'].fillna('neutral', inplace = True)
    df['impact'] = pd.Categorical(df['impact'], categories = ['LoF', 'neutral', 'GoF'], ordered = True)

    df['zscore'] = df['delta_score']/std
    df['zscore_bins'] = pd.cut(df['zscore'], bins = (-np.inf, -2,-1,0,1,2, np.inf)) 

    impact_per_sig = pd.pivot_table(df, index = 'impact', columns = 'INFO/CLNSIG',aggfunc = 'size').T
    impact_per_sig.div(impact_per_sig.sum(axis = 0)).T.plot.bar(
        cmap = 'coolwarm', figsize = (3,3))
    plt.ylabel('Fraction of variants')
    plt.legend(bbox_to_anchor = (1.5,1))
    sns.despine()
    plt.savefig(indir / 'output/variants' / 'clinvar_analysis' / f'{exp}.impact_CLNSIG.pdf')

    # parse by disease
    df['CLNDN_parsed'] = df['INFO/CLNDN'].str.split('|').tolist()
    df_exploded = df.explode('CLNDN_parsed')
    disease_impact_count = pd.pivot_table(df_exploded, index = 'CLNDN_parsed', 
                                        columns = 'INFO/CLNSIG', aggfunc = 'size').fillna(0)
    disease_impact_count.loc[disease_impact_count.sum(axis = 1).sort_values(ascending = False).index
        ].iloc[:10,:].plot.barh(stacked = True, figsize = (3,3), cmap = 'coolwarm')
    plt.xlabel('# variants overlapping \n binding site')
    plt.legend(title = 'CLNSIG', bbox_to_anchor = (1, 0.5))
    sns.despine()
    plt.savefig(indir / 'output/variants' / 'clinvar_analysis' / f'{exp}.variants_overlapping_binding.pdf')

    disease_impact_count = pd.pivot_table(df_exploded.loc[df['INFO/CLNSIG']!='Benign/Likely_benign'], index = 'CLNDN_parsed', 
                                      columns = 'zscore_bins', aggfunc = 'size').fillna(0)
    disease_impact_count.loc[disease_impact_count.sum(axis = 1).sort_values(ascending = False).index
    ].iloc[:10,:].plot.barh(stacked = True, figsize = (3,3), cmap = 'coolwarm')
    plt.xlabel('# non-benign variants')
    plt.legend(title = 'delta score bins')
    sns.despine()
    plt.savefig(indir / 'output/variants' / 'clinvar_analysis' / f'{exp}.non_benign_impact_cnt.pdf')

    disease_impact_count = pd.pivot_table(df_exploded.loc[df['INFO/CLNSIG']!='Benign/Likely_benign'], index = 'CLNDN_parsed', 
                                      columns = 'zscore_bins', aggfunc = 'size').fillna(0)
    disease_impact_count.div(disease_impact_count.sum(axis = 1
                                                    ), axis = 0).loc[disease_impact_count.sum(axis = 1).sort_values(ascending = False).index
        ].iloc[:10,:].sort_values(by = pd.Interval(-np.inf, -2, closed='right')).plot.barh(stacked = True, figsize = (3,3), cmap = 'coolwarm')
    sns.despine()
    plt.legend(bbox_to_anchor = (1.5,1))
    plt.xlabel('Fraction non-benign variants')
    plt.savefig(indir / 'output/variants' / 'clinvar_analysis' / f'{exp}.non_benign_impact_distribution.pdf')

    # # test for singificance for constrained genes
    gnomAD_test['is_ClinVar']=gnomAD_test.index.isin(df['gene_name'])
    gnomAD_test['sig_against_0'] = gnomAD_test['FDR_against_0']<0.2
    gnomAD_test['sig_against_mean'] = gnomAD_test['FDR_against_mean']<0.2

    cont = pd.pivot_table(gnomAD_test, index = 'is_ClinVar', columns = 'sig_against_0', aggfunc = 'size').fillna(0)
    pv, odds = testing(cont)
    with open(indir / 'output/variants' / 'clinvar_analysis' / f'{exp}.constrain_disease_test', 'w') as f:
        f.write(f'{exp},{pv},{odds}\n')

    # prioritize disease variants in constrained genes
    f, ax = plt.subplots(1,2, sharex = True, figsize = (12,3))
    lof_count = df.loc[(df['gene_name'].isin(gnomAD_test.loc[gnomAD_test['sig_against_0']].index))
    & (df['impact']=='LoF')
    ].explode('CLNDN_parsed')['CLNDN_parsed'].value_counts()
    if not lof_count.empty:
        lof_count.plot.barh(
            stacked = True, by = 'z_score_bins',
            ax = ax[0],
            color = 'royalblue'
        )
    ax[0].set_xlabel('# LoF variants in constrained genes')

    
    gof_count = df.loc[(df['gene_name'].isin(gnomAD_test.loc[gnomAD_test['sig_against_0']].index))
    & (df['impact']=='GoF')
    ].explode('CLNDN_parsed')['CLNDN_parsed'].value_counts()
    if not gof_count.empty:
        gof_count.plot.barh(
            stacked = True, by = 'z_score_bins',
            ax = ax[1],
            color = 'tomato'
        )
    ax[1].set_xlabel('# GoF variants in constrained genes')
    plt.tight_layout()
    sns.despine()

    disease_count = pd.concat([lof_count, gof_count], axis = 1).fillna(0)
    disease_count.columns = ['LoF', 'GoF']
    disease_count.to_csv(indir/'output/variants' / 'clinvar_analysis' / f'{exp}.constrain_disease_count.csv')


    f, ax = plt.subplots(1,2, sharex = True, figsize = (12,3))
    lof_count = df.loc[(df['gene_name'].isin(gnomAD_test.loc[gnomAD_test['sig_against_0']].index))
    & (df['impact']=='LoF')
    ]['gene_name'].value_counts()
    if not lof_count.empty:
        lof_count.plot.barh(
            stacked = True, by = 'z_score_bins',
            ax = ax[0],
            color = 'royalblue'
        )
    ax[0].set_xlabel('# LoF variants in constrained genes')

    gof_count = df.loc[(df['gene_name'].isin(gnomAD_test.loc[gnomAD_test['sig_against_0']].index))
    & (df['impact']=='GoF')
    ]['gene_name'].value_counts()
    if not gof_count.empty:
        gof_count.plot.barh(
            stacked = True, by = 'z_score_bins',
            ax = ax[1],
            color = 'tomato'
        )
    ax[1].set_xlabel('# GoF variants in constrained genes')
    plt.tight_layout()
    sns.despine()

    gene_count = pd.concat([lof_count, gof_count], axis = 1).fillna(0)
    gene_count.columns = ['LoF', 'GoF']
    disease_count.to_csv(indir/'output/variants' / 'clinvar_analysis' / f'{exp}.constrain_disease_gene_count.csv')

    # save file
    df.to_csv(indir/'output/variants' / 'clinvar' / f'{exp}.variants_scores.csv')
        

