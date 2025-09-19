import pandas as pd
from pathlib import Path
import sys
from pybedtools import BedTool
from scipy.stats import fisher_exact,chisquare
import numpy as np

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
    finemap = pd.read_csv(finemap_annotation, sep = '\t')
    tested = list((indir / 'output/tested_windows').glob(f'{exp}*'))
    annotation = pd.read_csv(sys.argv[3],
                            sep = '\t')
    outf = sys.argv[4]

    # additional files required
    top1 = BedTool('/tscc/nfs/home/hsher/ps-yeolab5/ukbb_dr_score/top1.bed')
    top5 = BedTool('/tscc/nfs/home/hsher/ps-yeolab5/ukbb_dr_score/top5.bed')

    # find all tested windows
    all_tested = set(annotation['name'])
    for f in tested:
        df = pd.read_csv(f, sep = '\t')
        all_tested = all_tested.intersection(set(df['name']))
    
    # gather information
    all_info = annotation.loc[annotation['name'].isin(all_tested)]
    all_info['is_binding'] = all_info['name'].isin(finemap['name.1'])

    # find top DR
    tested_bed = BedTool.from_dataframe(all_info[['chrom', 'start', 'end', 'name']])
    all_info['is_top1_DR_percentile'] = all_info['name'].isin(tested_bed.intersect(top1, f=1).to_dataframe()['name'].tolist())
    all_info['is_top5_DR_percentile'] = all_info['name'].isin(tested_bed.intersect(top5, f=1).to_dataframe()['name'].tolist())

    # test for enrichment of RBP binding site in topDR
    stats = []
    for col in ['is_top5_DR_percentile','is_top1_DR_percentile']:
        for name, group in all_info.groupby(by = 'feature_type_top'):
            pv = pd.pivot_table(group, index = 'is_binding', columns = col, aggfunc = 'size').fillna(0)
            pval, odds_ratio = testing(pv)
            stats.append([pval, odds_ratio, name, col])
    stats = pd.DataFrame(stats, columns = ['p-value', 'odds ratio', 'feature_type_top', 'subset'])

    stats.to_csv(outf)