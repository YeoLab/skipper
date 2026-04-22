import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from scipy.stats import chisquare
from scipy.stats import binom
def test_positive_selection(
    df,
    rate_index = ('neutral', '.')
    ):
    stat = []
    p_null = df.loc[rate_index, 'COSMIC_SAMPLE_MUTATED']/df.loc[rate_index, 'COSMIC_SAMPLE_TESTED']
    
    for index, row in df.iterrows():
        if index != rate_index:
            # probability of getting something
            pvalue = 1-binom.cdf(n = row['COSMIC_SAMPLE_TESTED'], k = row['COSMIC_SAMPLE_MUTATED'],
                                 p = p_null)
            fc =  (row['COSMIC_SAMPLE_MUTATED']/row['COSMIC_SAMPLE_TESTED'])/p_null

            stat.append(list(index)+[pvalue, fc])

    stat = pd.DataFrame(stat, columns = df.index.names + ['pvalue', 'fc'])
    return stat

def test_association(df,
    categories = ['LoF', 'GoF'], 
    groups = ['1','2'],
    null_group = '.'
                    ):
    '''
    impact_per_sig: categories by groups
    '''

    # get frequency of categories in group
    freq = df.div(df.sum(axis = 0), axis = 1)
    
                     
    stats = []
    null_frequency = freq[null_group]
    for category in categories: # LoF or GoF
        # s = maf_vs_impact.sum()
        for index in groups: # Tier 1,2
            try:
                row = df[index]
                total_cnt = row.sum()
                # expected count from tier '
                p = null_frequency.loc[category] # 
                
                f = np.array([p, 1-p])
                exp = total_cnt*f
        
                # observed
                cnt = row.loc[category]
                
                
                chi, p = chisquare(np.array([cnt, total_cnt-cnt]), exp)
                odds_ratio = (cnt/(total_cnt-cnt))/(p/(1-p))
                stats.append([index, category, p, odds_ratio])
            except Exception as e:
                print(e, index, category)
            
    stats = pd.DataFrame(stats, columns = ['group', 'category', 'pvalue', 'odds_ratio'])
    return stats

def sigmoid(z):
    return 1/(1 + np.exp(-z))

if __name__ == '__main__':
    indir = Path(sys.argv[1])
    exp = Path(sys.argv[2])
    mutation_census_file = '/tscc/nfs/home/hsher/ps-yeolab5/cosmic_data/CancerMutationCensus_AllData_v98_GRCh38.tsv.gz'
    gene_census_file = '/tscc/nfs/home/hsher/ps-yeolab5/cosmic_data/CancerGeneCensus.tsv'

    # load cancer gene annotations
    cancer_mut_census = pd.read_csv(mutation_census_file,
                                    sep = '\t')
    cancer_gene_census = pd.read_csv(gene_census_file,
                                    sep = '\t')
    cancer_gene_census['Role'] = cancer_gene_census['Role in Cancer'].str.split(', ')
    cancer_gene_census_exploded=cancer_gene_census.explode('Role')

    # calculate binding score
    ref = pd.read_csv(indir / 'output/variants' / 'cosmic_noncoding' / f'{exp}.ref.score.txt', sep = '\t', names = ['ID', 'score'])
    ref.rename({'score': 'score_ref'}, axis = 1, inplace = True)
    ref.drop('ID', axis = 1, inplace = True)
    alt = pd.read_csv(indir / 'output/variants' / 'cosmic_noncoding' / f'{exp}.alt.score.txt', sep = '\t', names = ['ID', 'score'])
    alt.rename({'score': 'score_alt'}, axis = 1, inplace = True)
    alt.drop('ID', axis = 1, inplace = True)
    annot = pd.read_csv(indir / 'output/variants' / 'cosmic_noncoding' / f'{exp}.csv', index_col = 0)

    df = pd.concat([ref,alt, annot], axis = 1)

    df.rename({'5': 'INFO/SAMPLE_COUNT', '6': 'INFO/TIER'}
            , axis = 1, inplace = True)

    df['score_alt'] = sigmoid(df['score_alt'])
    df['score_ref'] = sigmoid(df['score_ref'])
    df['delta_score']=df['score_alt']-df['score_ref']

    # annotate
    finemap_annotation = indir / 'output/finemapping/mapped_sites/' / f'{exp}.finemapped_windows.annotated.tsv'
    site_annotation_uniq = pd.read_csv(finemap_annotation, sep = '\t')


    # now annotate, all variants should be derived from finemapped windows
    assert df['name'].isin(site_annotation_uniq['name']).all()
    df['feature_type_top'] = df['name'].map(site_annotation_uniq.set_index('name')['feature_type_top'])
    df['transcript_type_top'] = df['name'].map(site_annotation_uniq.set_index('name')['transcript_type_top'])
    df['gene_name']=df['name'].map(site_annotation_uniq.set_index('name')['gene_name'])
    df['strand'] = df['name'].map(site_annotation_uniq.set_index('name')['strand'])

    # annotate variant type
    df.loc[(df['REF'].str.len()==1)&(df['ALT'].str.len()==1), 'TYPE']='SNV'
    df.loc[(df['REF'].str.len()!=1)|(df['ALT'].str.len()!=1), 'TYPE']='INDEL'
    df = df.merge(cancer_mut_census, left_on = 'ID', right_on = 'GENOMIC_MUTATION_ID')

    # annotate impact
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

    # annotate gene role
    gene2role = pd.pivot_table(cancer_gene_census_exploded,
               index = 'Gene Symbol',
               columns = ['Role'], aggfunc = 'size').fillna(0).applymap(lambda x: True if x==1 else False)
    df = df.merge(gene2role, left_on = 'gene_name', right_index = True, how = 'left')
    df[['oncogene', 'TSG']] = df[['oncogene', 'TSG']].fillna(False)

    # test association RBP impact with gene tier
    impact_per_sig = pd.pivot_table(df, index = 'impact', columns = 'INFO/TIER',aggfunc = 'size')
    impact_per_sig.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.impact_vs_tier.csv')
    
    stats_tier = test_association(impact_per_sig)
    stats_tier.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.impact_vs_tier_stat.csv')

    # test association RBP impact with gene role
    role_to_variants = df.merge(cancer_gene_census_exploded[['Gene Symbol', 'Role']],
         left_on = 'gene_name',
         right_on = 'Gene Symbol',
            how = 'left')
    role_to_variants['Role'].fillna('No established role', inplace = True)
    impact_per_role = pd.pivot_table(role_to_variants, index = 'impact', columns = 'Role',aggfunc = 'size')
    impact_per_role.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.impact_vs_role.csv')

    stats_role = test_association(impact_per_role, groups = ['TSG', 'fusion', 'oncogene'],
                             null_group = 'No established role')
    stats_role.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.impact_vs_role_stat.csv')

    # test positive selection in various groups
    saf_by_impact = df.groupby(by = ['impact', 'INFO/TIER'])['COSMIC_SAMPLE_MUTATED', 'COSMIC_SAMPLE_TESTED'].sum()
    saf_by_impact.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.SAF_by_tier.csv')
    try:
        stat_positive_selection_by_tier = test_positive_selection(saf_by_impact,
                                                                rate_index = ('neutral', '2'))

        
    except Exception as e:
        stat_positive_selection_by_tier = pd.DataFrame()
    
    stat_positive_selection_by_tier.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.pos_selection_tier.csv')
    
    try:
        saf_by_oncogene = df.groupby(by = ['impact', 'oncogene'])['COSMIC_SAMPLE_MUTATED', 'COSMIC_SAMPLE_TESTED'].sum()
        saf_by_oncogene.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.SAF_by_oncogene.csv')
        stat_positive_selection_by_oncogene = test_positive_selection(saf_by_oncogene, rate_index = ('neutral', True))
        
    except Exception as e:
        stat_positive_selection_by_oncogene = pd.DataFrame()
    stat_positive_selection_by_oncogene.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.pos_selection_oncogene.csv')

    try:
        saf_by_tsg = df.groupby(by = ['impact', 'TSG'])['COSMIC_SAMPLE_MUTATED', 'COSMIC_SAMPLE_TESTED'].sum()
        saf_by_tsg.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.SAF_by_oncogene.csv')
        stat_positive_selection_by_tsg = test_positive_selection(saf_by_tsg, rate_index = ('neutral', True))
    except Exception as e:
        stat_positive_selection_by_tsg = pd.DataFrame()
    stat_positive_selection_by_tsg.to_csv(indir/'output/variants' / 'cosmic_analysis' / f'{exp}.pos_selection_tsg.csv')

    df.to_csv((indir/'output/variants' / 'cosmic_noncoding' / f'{exp}.variants_scores.csv'))
