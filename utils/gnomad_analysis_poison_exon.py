import pandas as pd
import sys
import numpy as np
from scipy.stats import mannwhitneyu
import math
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

def make_MAF_bins(gnomad):
    gnomad['MAF']=gnomad['INFO/AC']/gnomad['INFO/AN']
    gnomad['MAF_bin'] = pd.cut(gnomad['MAF'], bins = [0, 1e-3, 1e-2, 1])
    
    
    bins = gnomad['MAF_bin'].cat.categories.tolist()
    bins = ['singleton']+bins
    gnomad['MAF_bin']=pd.Categorical(gnomad['MAF_bin'], categories = bins, ordered = True)
    gnomad.loc[gnomad['INFO/AC']==1, 'MAF_bin']='singleton'
    gnomad['MAF_bin_rank']=gnomad['MAF_bin'].rank(method = 'dense')
    gnomad['is_singleton']=gnomad['MAF_bin'].eq('singleton')
    return gnomad

if __name__ == '__main__':
    all_scores = pd.read_csv(sys.argv[1])
    gnomad = pd.read_csv(sys.argv[2],na_values = '.')
    out_prefix = sys.argv[3]


    std = gnomad.loc[(gnomad['INFO/AC']/gnomad['INFO/AN'])>0.01,'delta_score'].std()
    if math.isnan(std):
        std = gnomad['delta_score'].std()

    # poison_exon_annotation: 
    annotation = pd.read_csv('/tscc/nfs/home/hsher/bin/poison_exon_variant_analysis/data/Felker2023SuppelementaryTable1_hg38_PE_cassettes_and_elements_annotated.tsv',
                        sep = '\t')
    annotation['is_SFARI']=(~annotation['SFARI_genetic_category'].isnull())
    annotation['is_MIM']=(~annotation['MIM_phenotypes'].isnull())

    # gnomAD constraints
    gnomad_constrain = pd.read_csv('/tscc/projects/ps-yeolab5/hsher/gnomAD/gnomad.v4.0.constraint_metrics.tsv', sep = '\t')
    gnomad_constrain = gnomad_constrain.loc[(gnomad_constrain['mane_select'])&(gnomad_constrain['transcript'].str.contains('ENST'))]
    gnomad_constrain['lof.pLI bins']=pd.cut(gnomad_constrain['lof.pLI'], bins = [0,0.1,0.9,1],
                            labels = ['LoF tolerant', 'Ambiguous', 'LoF intolerant']
        
        )
    annotation['lof.pLI bins']=annotation['gene_symbol'].map(gnomad_constrain.dropna().set_index('gene')['lof.pLI bins'])

    # annotate impact
    all_scores['delta_score_bins'] = pd.cut(all_scores['delta_score'], bins = [-np.inf,-2*std, 2*std, np.inf],
                                )
    bin2name = dict(zip(all_scores['delta_score_bins'].cat.categories,
                        ['Weakening', 'No change', 'Strengthening']
                    )
                )
    reference_category = all_scores['delta_score_bins'].cat.categories[1]

    # annotate gene
    all_scores = all_scores.merge(annotation, left_on = 'name',
                 right_on = 'PE_identifier')
    
    # annotate MAF
    make_MAF_bins(all_scores)

    # load o/e reference
    reference_df= pd.read_csv('/tscc/nfs/home/hsher/ps-yeolab5/gnomAD_reference_set/v4/vep.csv',
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

    # scale to intergenic o/e
    scaling_rate = reference_df.loc['intergenic_variant']['o/e']

    # scale MAPS
    
    def sigmoid(x):
        return 1 / (1 + np.exp(-x))
    singleton_coef = pd.read_csv('/tscc/projects/ps-yeolab3/hsher/ENCODE/3_CLIP_ML_RBPNet/singleton_adjust.coef.csv', index_col = 0)['0']
    all_scores['expected_singleton']=sigmoid(
    all_scores['INFO/MR']*singleton_coef['MR']+singleton_coef['Intercept'])

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

    def add_oe_reference(ax, xmin, xmax, ylim, show_reference_label):
        for index, row in reference_df.loc[selected].iterrows():
            if row['o/e'] > ylim:
                ax.hlines(y = row['o/e'], xmin = -0.5, xmax = xmax, 
                        color = 'lightgrey', linestyle = 'dotted', label = index)
                if show_reference_label:
                    ax.text(x = xmax, y = row['o/e'], s = index)

    
    f, ax = plt.subplots(figsize = (4,4))
    bydelta_ind = groupby_bootstrap(all_scores, by = ['delta_score_bins'], n=10)

    x = 'delta_score_bins'
    y = 'o/e'
    global_oe = oe_test_for_negative_selection(bydelta_ind, x, y)
    global_oe.to_csv(f'{out_prefix}.global_oe.csv')
    
    sub_df = all_scores.loc[~all_scores['INFO/AN'].isnull()]
    bydelta_maps = groupby_bootstrap_maps(sub_df, by = ['delta_score_bins'])
    x = 'delta_score_bins'
    y = 'MAPS'
    global_maps = oe_test_for_negative_selection(bydelta_maps, x, y)
    global_maps.to_csv(f'{out_prefix}.global_maps.csv')

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
    
    by = 'is_MIM'
    bydelta, feature_results = groupby_test_results(all_scores, by, y = 'o/e')
    feature_results.to_csv(f'{out_prefix}.MIM_oe.csv')

    by = 'is_SFARI'
    bydelta, feature_results = groupby_test_results(all_scores, by, y = 'o/e')
    feature_results.to_csv(f'{out_prefix}.SFARI_oe.csv')

    by = 'lof.pLI bins'
    bydelta, feature_results = groupby_test_results(all_scores, by, y = 'o/e')
    feature_results.to_csv(f'{out_prefix}.LoF_oe.csv')




