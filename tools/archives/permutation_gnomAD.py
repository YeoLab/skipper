import pandas as pd
import statsmodels.formula.api as smf
import concurrent.futures
import random
import time
import multiprocessing
import sys
import matplotlib.pyplot as plt
from scipy.stats import norm, uniform
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.api as sm


def lm(subset_df):
    ''' test association between delta score and MAF'''
    results = smf.ols('delta_score ~ MAF_bin_rank', data=subset_df).fit()
    return results.params['MAF_bin_rank']

def groupby_lm(df, by = 'gene_name', max_workers = 8):
    ''' for each "by", test association '''
    keys = [name for name, group in df.groupby(by = by)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        result = list(executor.map(lm, [group for name, group in df.groupby(by = by)]))
    
    return dict(zip(keys, result))

def permute(group, value_to_permute = 'delta_score'):
    val = group[value_to_permute].tolist()
    random.shuffle(val)
    group[value_to_permute] = val
    return group

def groupby_permute(df, groupby_col = 'MAF_bin_rank', max_workers = 8, by = 'gene_name'):
    ''' permute to get null distribution per "by" column '''

    print('start permuting')
    grouped = df.groupby(groupby_col)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        resultlist = list(executor.map(permute, [group for name, group in grouped]))
    permuted_df = pd.concat(resultlist, axis = 0)
    print('end permuting, starting lm')
    
    result = groupby_lm(permuted_df, max_workers = max_workers, by = by)
    print('done')
    del permuted_df
    return result

if __name__ == '__main__':

    n_cpu = multiprocessing.cpu_count()
    print(f'Number of CPU available: {n_cpu}')
    col_to_test = sys.argv[2]
    n_iter=100
    
    # input must have 'MAF_bin_rank', 'delta_score' and 'TYPE'
    df = pd.read_csv(sys.argv[1],
                    sep = '\t', comment = '#',
                    usecols=['MAF_bin_rank', 'delta_score', col_to_test, 'TYPE'])
    
    out_prefix = sys.argv[3]
    
    # only test for SNV
    df = df.loc[df['TYPE']=='SNV']
    df.drop('TYPE', inplace = True, axis = 1)

    assert not df['delta_score'].isnull().any()

    # generate null distribution
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        result = executor.map(groupby_permute, [df]*n_iter)
    
    null = pd.DataFrame(list(result))
    null.to_csv(f'{out_prefix}.null.csv')

    # do the test
    mean_null, std_null = null.mean(axis = 0),null.std(axis = 0)
    stat = pd.concat([mean_null, std_null], axis = 1)
    stat.columns = ['mean', 'std']

    # get slope
    slope = groupby_lm(df)
    stat['slope'] = stat.index.map(slope)
    
    # test against slope>0
    stat['pvalue_against_0'] = stat.apply(lambda row: 
                            1-norm(loc = 0, scale = row['std']).cdf(row['slope']),
                         axis = 1)
    _, stat['FDR_against_0'] = fdrcorrection(stat['pvalue_against_0'], alpha = 0.2)

    # test against slope>avg (transcriptome wide mean)
    stat['pvalue_against_mean'] = stat.apply(lambda row: 
                            1-norm(loc = row['mean'], scale = row['std']).cdf(row['slope']),
                         axis = 1)
    _, stat['FDR_against_mean'] = fdrcorrection(stat['pvalue_against_mean'], alpha = 0.2)

    # plott qq plot
    f, ax = plt.subplots(1,2)
    p = sm.qqplot(stat['pvalue_against_mean'], line ='45',
             dist=uniform, ax = ax[0])
    ax[0].set_title('H1: slope > transcriptome-wide)')

    p = sm.qqplot(stat['pvalue_against_0'], line ='45',
             dist=uniform, ax = ax[1])
    ax[1].set_title('H1: slope > 0)')
    plt.savefig(f'{out_prefix}.qqplot.pdf')

    # savefile
    stat.to_csv(f'{out_prefix}.test_stats.csv')
