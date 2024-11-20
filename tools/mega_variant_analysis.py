import pandas as pd
from pathlib import Path
import pybedtools
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import fisher_exact,chisquare,mannwhitneyu
import statsmodels.api as sm
import math
from statsmodels.stats.multitest import fdrcorrection
import warnings

def make_MAF_bins(gnomad):
    """ create MAF bins"""
    gnomad['MAF']=gnomad['INFO/AC']/gnomad['INFO/AN']
    gnomad['MAF_bin'] = pd.cut(gnomad['MAF'], bins = [0, 1e-3, 1e-2, 1])
    
    
    bins = gnomad['MAF_bin'].cat.categories.tolist()
    bins = ['singleton']+bins
    gnomad['MAF_bin']=pd.Categorical(gnomad['MAF_bin'], categories = bins, ordered = True)
    gnomad.loc[gnomad['INFO/AC']==1, 'MAF_bin']='singleton'
    gnomad['MAF_bin_rank']=gnomad['MAF_bin'].rank(method = 'dense')

    return gnomad

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

def test_subset(df, col = 'feature_type_top'):
    ''' test enrichment in subgroups based on columns'''
    all_stats = []
    
    for name, group in df.groupby(by = col):
        if group.shape[0]>100:
            maf_vs_impact = pd.pivot_table(group,                
                        index = 'MAF_bin', 
                       columns = 'impact', 
                       aggfunc = 'size')
            maf_vs_impact.div(maf_vs_impact.sum(axis = 0)).T.plot.bar(
                cmap = 'coolwarm', figsize = (3,3))
            plt.ylabel('Fraction of variants')
            plt.legend(bbox_to_anchor = (1.5,1))
            plt.title(name)
            sns.despine()
            plt.show()
            stats_local = test_enrichment_outliers(maf_vs_impact)
            stats_local[col]=name
            all_stats.append(stats_local)
    all_stats = pd.concat(all_stats,axis = 0)
    all_stats['sig'], all_stats['FDR']= fdrcorrection(all_stats['pvalue'], alpha = 0.2)
    return all_stats
def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"
def sigmoid(x):
  return 1 / (1 + np.exp(-x))
if __name__ == '__main__':
    indir = Path(sys.argv[1])
    exp = sys.argv[2]
    outdir = indir / 'output/variant_analysis'

    # read gnomAD and roulette
    df = pd.read_csv(indir / f'output/variants/gnomAD_roulette/{exp}.total.csv',
                na_values = '.')
    annotation = pd.read_csv(indir / f'output/finemapping/mapped_sites/{exp}.finemapped_windows.annotated.tsv', sep = '\t')

    # define std from commmon variants
    df = make_MAF_bins(df)
    std = df.loc[df['MAF']>0.01, 'delta_score'].std()
    if math.isnan(std):
        warnings.warn('No common variants observed in binding site. Using all variant to estimate std')
        std = df['delta_score'].std()
    df.loc[df['delta_score']<-2*std, 'impact'] = 'LoB'
    df.loc[df['delta_score']>2*std, 'impact'] = 'GoB'
    df['impact'].fillna('neutral', inplace = True)

    # annotate 
    df['feature_type_top'] = df['name'].map(annotation.set_index('name')['feature_type_top'])
    df['transcript_type_top'] = df['name'].map(annotation.set_index('name')['transcript_type_top'])
    df['gene_name'] = df['name'].map(annotation.set_index('name')['gene_name'])

    gnomad_constrain = pd.read_csv('/tscc/projects/ps-yeolab5/hsher/gnomAD/gnomad.v4.0.constraint_metrics.tsv', sep = '\t')
    gnomad_constrain = gnomad_constrain.loc[(gnomad_constrain['mane_select'])&(gnomad_constrain['transcript'].str.contains('ENST'))]
    gnomad_constrain['lof.pLI bins']=pd.cut(gnomad_constrain['lof.pLI'], bins = [0,0.1,0.9,1],
                            labels = ['LoF tolerant', 'Ambiguous', 'LoF intolerant']
        
        )
    df['lof.pLI bins']=df['gene_name'].map(gnomad_constrain.dropna().set_index('gene')['lof.pLI bins'])
    df['mis_pphen.oe']=df['gene_name'].map(gnomad_constrain.dropna().set_index('gene')['mis_pphen.oe'])
    df['syn.oe']=df['gene_name'].map(gnomad_constrain.dropna().set_index('gene')['syn.oe'])

    # test enrichment of rare variants
    maf_vs_impact = pd.pivot_table(df,                
            index = 'MAF_bin', 
           columns = 'impact', 
           aggfunc = 'size')
    test_enrichment_outliers(maf_vs_impact).to_csv(outdir / f'{exp}.global_spectrum_enrichment.csv')

    for col in ['feature_type_top', 'transcript_type_top']:
        test_subset(df, col = col).to_csv(outdir / f'{exp}.{col}_spectrum_enrichment.csv')
    
    # load mutation rate expected singleton ratio for MAPs
    # singleton_model = sm.load('/tscc/nfs/home/hsher/projects/ENCODE/3_CLIP_ML_RBPNet/singleton_adjust.pickle') ## TODO: this has to be outside
    singleton_coef = pd.read_csv('/tscc/nfs/home/hsher/projects/ENCODE/3_CLIP_ML_RBPNet/singleton_adjust.coef.csv', index_col = 0)['0']
    df['expected_singleton']=sigmoid(df['INFO/MR']*singleton_coef['MR']+singleton_coef['Intercept'])
    df['is_singleton']=df['MAF_bin'].eq('singleton')
    

    # o/e ratio reference
    reference_df= pd.read_csv('/tscc/nfs/home/hsher/ps-yeolab5/gnomAD_reference_set/v4/vep.csv', ## TODO: this has to be outside
                        index_col = 0
                       )
    selected = ['3_prime_UTR_variant',
    '5_prime_UTR_variant',
    'intron_variant',
    'missense_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'synonymous_variant']

    reference_df.loc[selected].sort_values(by = 'o/e').reset_index().plot.scatter(
        x = 'Consequence_splitted', y = 'o/e', color = 'grey', figsize = (3,3))
    plt.xticks(rotation = 90)
    sns.despine()

    # scale to intergenic o/e
    scaling_rate = reference_df.loc['intergenic_variant']['o/e']

    # delta_score
    df['delta_score_bins'] = pd.cut(df['delta_score'], bins = [-np.inf,-2*std, 2*std, np.inf],)
    df['delta_zscore'] = pd.cut(df['delta_score'], bins = [-np.inf,-2*std,-std, std ,2*std, np.inf],
                            labels = ['LoB', 'weak LoB', 'Neutral', 'weak GoB', 'GoB']
                               )
    print(df['delta_score_bins'].value_counts())
    bin2name = dict(zip(df['delta_score_bins'].cat.categories,
                    ['Weakening', 'No change', 'Strengthening']
                   )
               )
    reference_category = df['delta_score_bins'].cat.categories[1]
    sub_df = df.loc[~df['INFO/AN'].isnull()]
    df.to_csv(outdir / f'{exp}.annotated.csv.gz')

    # bootstraping function for o/e and MAPS
    def groupby_bootstrap(df, by = ['delta_score_bins'], n = 30, ratio = 0.8, scaling_rate = scaling_rate):
        ''' groupby bootstrap for o/e '''
        
        data = []
        for name, group in df.groupby(by = by):
            sample_size = group.shape[0]
            for i in range(n):
                n_sampled = int(sample_size*ratio)
                sampled = group.sample(n_sampled)
                sampled_oe = (sampled['MAF'].count()/sampled['INFO/MR'].sum())/scaling_rate
                if type(name)==tuple:
                    data.append(list(name)+[sampled_oe, n_sampled, sample_size])
                else:
                    data.append([name, sampled_oe, n_sampled, sample_size])
        data = pd.DataFrame(data, columns = by+['o/e','n', 'total'])
        return data
    def groupby_bootstrap_maps(df, by = ['delta_score_bins'], n = 30, ratio = 0.8):
        data = []
        for name, group in df.groupby(by = by):
            sample_size = group.shape[0]
            for i in range(n):
                n_sampled = int(sample_size*ratio)
                sampled = group.sample(n_sampled)
                n_singleton=sampled['is_singleton'].sum()
                expected_singleton=sampled['expected_singleton'].sum()
                maps = (n_singleton-expected_singleton)/n_sampled
                
                
                if type(name)==tuple:
                    data.append(list(name)+[maps, n_sampled, sample_size])
                else:
                    data.append([name, maps, n_sampled, sample_size])
        data = pd.DataFrame(data, columns = by+['MAPS','n', 'total'])
        return data
    import math
    def oe_test_for_negative_selection(bydelta_ind, x, y, reference = reference_category,
                                    n_variants = 30):
        n_category = bydelta_ind[x].unique().shape[0]
        if reference is None:
            #use the middle class as reference
            reference = bydelta_ind[x].unique()[math.floor(n_category/2)]

        if y == 'o/e':
            alternative = 'less'
        elif y == 'MAPS':
            alternative = 'greater'
        stat_results = []
        for i, (name, group) in enumerate(bydelta_ind.groupby(by = x)):
            
            if name != reference and group['total'].max()>n_variants: # require at least 30 variants to test
                try:
                    _, p = mannwhitneyu(group[y], bydelta_ind.loc[bydelta_ind[x]==reference, y],
                                    alternative = alternative) # test for negative selection only
                    median_diff = group[y].median()-bydelta_ind.loc[bydelta_ind[x]==reference, y].median()
                    stat_results.append([bin2name[name], p, median_diff, group['total'].max()])
                except Exception as e:
                    print(e)
        return pd.DataFrame(stat_results, columns = [x, 'p-value', 'median difference', 'n'])


    ### global MAPS ###
    global_maps = groupby_bootstrap_maps(sub_df)
    global_result = oe_test_for_negative_selection(global_maps, 'delta_score_bins', 'MAPS')
    global_maps.to_csv(outdir / f'{exp}.global_MAPS.csv')
    global_result.to_csv(outdir / f'{exp}.global_MAPS_stat.csv')

    ### global o/e ###
    global_oe = groupby_bootstrap(df, by = ['delta_score_bins'])
    global_oe_result  = oe_test_for_negative_selection(global_oe, 'delta_score_bins', 'o/e')
    global_oe.to_csv(outdir / f'{exp}.global_oe.csv')
    global_oe_result.to_csv(outdir / f'{exp}.global_oe_stat.csv')

    ### testing subset of variants ###
    def groupby_test_results(df, by, y = 'o/e', n = 30):
        if type(by) == str:
            by = [by]

        if y == 'o/e':
            bydelta = groupby_bootstrap(df, by = ['delta_score_bins']+by, n=n)
        elif y == 'MAPS':
            bydelta = groupby_bootstrap_maps(df, by = ['delta_score_bins']+by, n=n)
        stats = []
        
        for fname, fgroup in bydelta.groupby(by = by):
            test_result = oe_test_for_negative_selection(fgroup, x='delta_score_bins', y=y)
            if type(fname) == tuple:
                for i, col in enumerate(by):
                    test_result[col]=fname[i]
            else:
                test_result[by]=fname
            stats.append(test_result)
        results = pd.concat(stats, axis = 0)

        return bydelta, results
    
    ### feature/transcript types ###
    for col in ['feature_type_top', 'transcript_type_top']:
        for metric, metric_name in zip(['o/e', 'MAPS'], ['oe', 'MAPS']):
            if metric == 'MAPS':
                bootstrap_data, test_results = groupby_test_results(sub_df, col, y = metric)
            else:
                bootstrap_data, test_results = groupby_test_results(df, col, y = metric)
            bootstrap_data.to_csv(outdir / f'{exp}.{col}.{metric_name}.csv')
            test_results.to_csv(outdir / f'{exp}.{col}.{metric_name}_stat.csv')
    
    ### stratify by LoF bins ###
    col = 'feature_type_top'
    for metric, metric_name in zip(['o/e', 'MAPS'], ['oe', 'MAPS']):
        if metric == 'MAPS':
            bootstrap_data, test_results = groupby_test_results(sub_df, ['lof.pLI bins', col],
                                                                    y = metric)
        else:
            bootstrap_data, test_results = groupby_test_results(df, ['lof.pLI bins', col],
                                                                    y = metric)
        bootstrap_data.to_csv(outdir / f'{exp}.{col}.{metric_name}.lofbins.csv')
        test_results.to_csv(outdir / f'{exp}.{col}.{metric_name}_stat.lofbins.csv')
    
    ### ClinVar ###
    clinvar_vcf = pd.read_csv(indir / 'output/variants' / 'clinvar' / f'{exp}.vcf', index_col = 0,
                         names = ['CHROM', 'POS', 'Clinvar_ID', 'REF', 'ALT',
                                  'INFO/CLNDN','INFO/CLNVC','INFO/CLNSIG',
                                  'INFO/CLNDISDB','INFO/AF_ESP','INFO/AF_EXAC',
                                  'INFO/AF_TGP','INFO/ALLELEID'],
    sep = '\t'
                            )

    clinvar_df = df.merge(clinvar_vcf, left_on = ['CHROM', 'POS', 'REF', 'ALT'],
                    right_on = ['CHROM', 'POS', 'REF', 'ALT']
                    )
    clinvar_df.drop_duplicates(subset = ['CHROM', 'POS', 'REF', 'ALT'],
                    inplace = True)
    
    clinvar_vep = pd.read_csv(indir / 'output/variants' / 'clinvar' / f'{exp}.vep.tsv', sep = '\t',
                          comment = '#',
                        names = ['Uploaded_variation','Location','Allele',
                          'Gene','Feature',
                          'Feature_type','Consequence','cDNA_position','CDS_position',
                          'Protein_position','Amino_acids','Codons','Existing_variation',
                          'Extra'])
    #clinvar_dedup = clinvar_vep[['Uploaded_variation', 'Location', 'Allele']].drop_duplicates()
    clinvar_df['Consequence'] = clinvar_df['Clinvar_ID'].map(
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
    def is_coding(s):
        if 'missense' in s or 'stop' in s or 'start' in s or 'frame' in s:
            return True
        else:
            return False
    clinvar_df['is_coding']=clinvar_df['Consequence'].apply(lambda x: is_coding(x))
    clinvar_df['CLNDN_parsed'] = clinvar_df['INFO/CLNDN'].str.split('|').tolist()
    clinvar_df.to_csv(outdir / f'{exp}.clinvar_variants.csv')
    
    # perform counting
    clinvar_df_exploded = clinvar_df.explode('CLNDN_parsed')
    clinvar_df_exploded.to_csv(outdir / f'{exp}.clinvar_variants_exploded.csv')

    def plot_sorted(counts, ax, n=10, **kwargs):
        counts.loc[counts.sum(axis = 1).sort_values(ascending = False).index
        ].iloc[:n,:].plot.barh(stacked = True, **kwargs, ax = ax)
    
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
    f, ax = plt.subplots(2,1, figsize = (6,6))
    try:
        plot_sorted(impact_count[True], ax = ax[0], cmap = 'coolwarm', legend = False)
    except KeyError as e:
        print(e)
    ax[0].set_title('coding')
    try:
        plot_sorted(impact_count[False], ax = ax[1], cmap = 'coolwarm')
    except KeyError as e:
        print(e)
    ax[1].set_title('non-coding')
    ax[1].set_xlabel('Number Variants')
    sns.despine()
    plt.savefig(outdir / f'{exp}.clinvar_impact_counts.pdf', bbox_inches='tight')
    impact_count.to_csv(outdir / f'{exp}.clinvar_impact_counts.csv')
