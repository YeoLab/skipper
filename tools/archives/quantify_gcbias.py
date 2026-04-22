import sys
import pandas as pd
from pathlib import Path
import numpy as np

def get_gc_odds_ratio(gc_bin_sum, samp, ip_rep, in_rep):
    gc_odds = (gc_bin_sum[f'{samp}_IP_{ip_rep}']/gc_bin_sum[f'{samp}_IN_{in_rep}']
          )/(total[f'{samp}_IP_{ip_rep}']/total[f'{samp}_IN_{in_rep}'])

    return gc_odds
def gc_fraction(gc_bin_sum, samp, ip_rep, in_rep):
    gc_frac = gc_bin_sum[f'{samp}_IP_{ip_rep}']/(gc_bin_sum[f'{samp}_IN_{in_rep}']+gc_bin_sum[f'{samp}_IP_{ip_rep}'])
    
    return gc_frac
def get_gc_odds_ratio_across_rep(gc_bin_sum, samp, lib = 'IP'):
    gc_odds = (gc_bin_sum[f'{samp}_{lib}_1']/gc_bin_sum[f'{samp}_{lib}_2']
          )/(total[f'{samp}_{lib}_1']/total[f'{samp}_{lib}_2'])

    return gc_odds

if __name__ == '__main__':
    counts_f = Path(sys.argv[1])
    samp = sys.argv[2]
    outf = sys.argv[3]
    counts = pd.read_csv(counts_f, sep = '\t')
    counts['gc_bin']=pd.qcut(counts['gc'], q = 10)

    bias_data = {}

    # calculate gcbias
    gc_bin_sum = counts.groupby(by = 'gc_bin')[counts.columns[7:-1]].sum()
    total = gc_bin_sum.sum(axis = 0)

    ip_reps = [1,2]
    in_reps = [1,1] if counts.columns.str.contains('_IN_').sum()==1 else [1,2]
    gc_odds_across_rep = []
    for ip_rep, in_rep in zip(ip_reps, in_reps):
        gc_odds = get_gc_odds_ratio(gc_bin_sum, samp, ip_rep, in_rep)
        gc_odds_across_rep.append(gc_odds)
    
    gc_bias = pd.concat(gc_odds_across_rep, axis = 1)
    gc_bias.columns = [f'Odds Ratio in IP_1', f'Odds Ratio in IP_2']
    bias_data['ip_in_l2or_diff'] = (np.log2(gc_bias['Odds Ratio in IP_1'])-np.log2(gc_bias['Odds Ratio in IP_2'])).abs().sum()

    ip_rep_bias = get_gc_odds_ratio_across_rep(gc_bin_sum, samp)
    bias_data['ip_l2odds_diff'] = np.log2(ip_rep_bias).abs().sum()

    if in_reps == [1,2]:
        in_rep_bias = get_gc_odds_ratio_across_rep(gc_bin_sum, samp, lib = 'IN')
        bias_data['in_l2odds_diff'] = np.log2(in_rep_bias).abs().sum()
    
    pd.Series(bias_data,name = samp).to_csv(outf)