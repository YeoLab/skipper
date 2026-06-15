import pandas as pd
from pathlib import Path
from pybedtools import BedTool
import sys
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams["image.cmap"] = "Dark2"
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Dark2.colors)
import seaborn as sns
from scipy.stats import norm
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection

def sigmoid(z):
    return 1/(1 + np.exp(-z))

def test_selection_in_subset(df, groupby = 'feature_type_top'):
    stats = []
    for name, group in df.groupby(by = groupby):
        results = smf.ols('delta_score ~ MAF_bin_rank', data=group).fit()
        stats.append([name, results.pvalues['MAF_bin_rank'], results.params['MAF_bin_rank'], group.shape[0]])
    stats = pd.DataFrame(stats, columns = [groupby,'pvalue', 'coef', '# SNP in binding site'])
    stats['sig'], stats['FDR']= fdrcorrection(stats['pvalue'], alpha = 0.2)

    return stats

if __name__ == '__main__':
    indir = Path(sys.argv[1])
    exp = sys.argv[2]
    finemap_annotation = indir / 'output/finemapping/mapped_sites/' / f'{exp}.finemapped_windows.annotated.tsv'

    # load all files
    ref = pd.read_csv(indir / 'output/variants' / 'gnomAD' / f'{exp}.ref.score.txt', sep = '\t', names = ['ID', 'score'])
    ref.rename({'score': 'score_ref'}, axis = 1, inplace = True)
    ref.drop('ID', axis = 1, inplace = True)
    alt = pd.read_csv(indir / 'output/variants' / 'gnomAD' / f'{exp}.alt.score.txt', sep = '\t', names = ['ID', 'score'])
    alt.rename({'score': 'score_alt'}, axis = 1, inplace = True)
    alt.drop('ID', axis = 1, inplace = True)
    annot = pd.read_csv(indir / 'output/variants' / 'gnomAD' / f'{exp}.csv', index_col = 0)
    site_annotation_uniq = pd.read_csv(finemap_annotation, sep = '\t')

    # assure they are the same length
    assert ref.shape[0] == alt.shape[0]
    assert ref.shape[0] == annot.shape[0]

    # calculate MAF and binning
    df = pd.concat([ref,alt, annot], axis = 1)
    df.rename({'5': 'INFO/AC', '6': 'INFO/AN'}, axis = 1, inplace = True)
    df = df.loc[df['INFO/AC']>0]
    # calculate scoring: sigmoid transformation into probability
    df['score_alt'] = sigmoid(df['score_alt'])
    df['score_ref'] = sigmoid(df['score_ref'])

    df['delta_score']=df['score_alt']-df['score_ref']
    df ['log_ratio'] = np.log(df['score_alt'])-np.log(df['score_ref'])

    df['MAF']=df['INFO/AC']/df['INFO/AN']
    df['MAF_bin'] = pd.cut(df['MAF'], bins = [0, 1e-3, 1e-2, 1])
    bins = df['MAF_bin'].cat.categories.tolist()
    bins = ['singleton']+bins
    df['MAF_bin']=pd.Categorical(df['MAF_bin'], categories = bins, ordered = True)
    df.loc[df['INFO/AC']==1, 'MAF_bin']='singleton'
    df['MAF_bin_rank']=df['MAF_bin'].rank(method = 'dense')

    

    # now annotate, all variants should be derived from finemapped windows
    assert df['name'].isin(site_annotation_uniq['name']).all()
    df['feature_type_top'] = df['name'].map(site_annotation_uniq.set_index('name')['feature_type_top'])
    df['transcript_type_top'] = df['name'].map(site_annotation_uniq.set_index('name')['transcript_type_top'])
    df['gene_name']=df['name'].map(site_annotation_uniq.set_index('name')['gene_name'])
    df['strand'] = df['name'].map(site_annotation_uniq.set_index('name')['strand'])

    # annotate variant type
    df.loc[(df['REF'].str.len()==1)&(df['ALT'].str.len()==1), 'TYPE']='SNV'
    df.loc[(df['REF'].str.len()!=1)|(df['ALT'].str.len()!=1), 'TYPE']='INDEL'

    # Summary statistics: Distribution of delta
    std=df['delta_score'].std()
    mean = df['delta_score'].mean()
    rv = norm(loc = 0, scale = std)
    
    f, ax = plt.subplots(2,1, sharex = True)
    
    x = np.linspace(df['delta_score'].min(),
                df['delta_score'].max(), 100)

    ax[0].plot(x, rv.pdf(x), 'r-', lw=2, label='Normal distribution\n(mean=0)')
    df['delta_score'].plot.hist(bins = 100, color = 'lightgrey', density = True, ax = ax[0])
    ax[0].vlines(x = std*2, ymin = 0, ymax = 1, color = 'black', 
            linestyle='dashed')
    ax[0].vlines(x = -std*2, ymin = 0, ymax = 1, color = 'black', 
            linestyle='dashed', label = f'2 STDEV(std={std:.2f})')
    ax[0].set_xlabel('delta score')
    ax[0].legend()

    sns.histplot(data = df, x = 'delta_score',
                     bins = 100, color = 'lightgrey', hue = 'TYPE',
            stat = 'frequency', ax = ax[1])
    plt.savefig(indir / 'output/variants' / 'gnomAD_analysis' / f'{exp}.delta_score_distribution.pdf')

    # annotate impact
    df.loc[df['delta_score']<-2*std, 'impact'] = 'LoF'
    df.loc[df['delta_score']>2*std, 'impact'] = 'GoF'
    df['impact'].fillna('neutral', inplace = True)

    # Summary statistics: MAF vs impact
    maf_vs_impact = pd.pivot_table(df,                
                index = 'MAF_bin', 
               columns = 'impact', 
               aggfunc = 'size')


    maf_vs_impact.div(maf_vs_impact.sum(axis = 0)).T.plot.bar(
        cmap = 'coolwarm', figsize = (3,3))
    plt.ylabel('Fraction of variants')
    plt.legend(bbox_to_anchor = (1.5,1))
    sns.despine()
    plt.savefig(indir / 'output/variants' / 'gnomAD_analysis' / f'{exp}.MAF_vs_impact.pdf')

    # test by subset
    feature_stat = test_selection_in_subset(df, groupby = 'feature_type_top')
    transcript_stat = test_selection_in_subset(df, groupby = 'transcript_type_top')

    feature_stat.to_csv(indir / 'output/variants' / 'gnomAD_analysis' / f'{exp}.feature_stat.tsv',
        sep = '\t')
    transcript_stat.to_csv(indir / 'output/variants' / 'gnomAD_analysis' / f'{exp}.transcript_stat.tsv',
        sep = '\t')
    
    # plot by subset
    f, ax = plt.subplots(1,2, figsize = (16,8), sharey = True)

    sns.boxplot(data = df,
                hue = 'MAF_bin',
                y = 'delta_score',
                x = 'feature_type_top',
                showfliers=False,
                ax = ax[0],
                palette="Blues"
                )
    ax[0].tick_params(axis='x', rotation=90)
    ax[0].get_legend().remove()

    sns.boxplot(data = df,
                hue = 'MAF_bin',
                y = 'delta_score',
                x = 'transcript_type_top',
                showfliers=False,
                ax = ax[1],
                palette="Blues"
                )
    ax[1].tick_params(axis='x', rotation=90)
    sns.move_legend(ax[1], "upper left", bbox_to_anchor=(1, 0.5), title='MAF')
    plt.tight_layout()
    sns.despine()
    plt.savefig(indir / 'output/variants' / 'gnomAD_analysis' / f'{exp}.subset_selection_boxplot.pdf')

    # save results
    outf = indir / 'output/variants' / 'gnomAD' / f'{exp}.variants_scores.tsv'

    with open(outf, 'w') as f:
        f.write(f'#std:{std}\n#mean: {mean}\n')
        df.to_csv(f, sep = '\t')
    
