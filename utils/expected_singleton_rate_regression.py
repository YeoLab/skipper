import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import statsmodels.api as sm
import statsmodels.formula.api as smf
def make_MAF_bins(gnomad):
    gnomad['MAF']=gnomad['INFO/AC']/gnomad['INFO/AN']
    gnomad['MAF_bin'] = pd.cut(gnomad['MAF'], bins = [0, 1e-3, 1e-2, 1],
                              labels = ['ultrarare',
                                        'rare',
                                        'common']
                                        
                              )
    
    
    bins = gnomad['MAF_bin'].cat.categories.tolist()
    bins = ['singleton']+bins
    gnomad['MAF_bin']=pd.Categorical(gnomad['MAF_bin'], categories = bins, ordered = True,
                                    )
    gnomad.loc[gnomad['INFO/AC']==1, 'MAF_bin']='singleton'
    gnomad['MAF_bin_rank']=gnomad['MAF_bin'].rank(method = 'dense')

    return gnomad
if __name__ == '__main__':
    roulette = pd.read_csv(sys.argv[1],
                       sep = '\t',
                       names = ['CHROM','POS', 'ID','REF','ALT', 'INFO/RBPSITE','INFO/MR', 'INFO/AC','INFO/AN'],
                       na_values = '.',
                      )
    vep = pd.read_csv(sys.argv[2],
                 sep = '\t',
                 names = ['Uploaded_variation','Location','Allele',
                          'Gene','Feature',
                          'Feature_type','Consequence','cDNA_position','CDS_position',
                          'Protein_position','Amino_acids','Codons','Existing_variation',
                          'Extra'])
    model_out = Path(sys.argv[3])

    vep = vep.loc[vep['Consequence']=='intergenic_variant']
    vep['POS']=vep['Location'].str.split(':', expand = True)[1].astype(int)
    merged_vep = vep.merge(roulette, right_on = ['POS', 'ALT'],
             left_on = ['POS', 'Allele'], how = 'left'
            )
    merged_vep = make_MAF_bins(merged_vep)
    merged_vep['mutation rate bin'] = pd.cut(merged_vep['INFO/MR'], bins = 100,
                                         duplicates = 'drop')
    count_per_mr_maf = pd.pivot_table(merged_vep, index = 'mutation rate bin', 
                                  columns = 'MAF_bin',
              aggfunc = 'size').dropna()
    fraction_per_mr_maf = count_per_mr_maf.div(count_per_mr_maf.sum(axis =1),
                                           axis = 0)
    types = fraction_per_mr_maf.columns
    fraction_per_mr_maf['MR']=[(i.left+i.right)/2 for i in fraction_per_mr_maf.index]
    fraction_per_mr_maf['counts']=merged_vep['mutation rate bin'].value_counts()
    
    fraction_per_mr_maf['log_counts'] = np.log10(fraction_per_mr_maf['counts'])
    count_per_mr_maf['MR']=[(i.left+i.right)/2 for i in count_per_mr_maf.index]
    for t in types:
        count_per_mr_maf[f'non_{t}']=count_per_mr_maf.sum(axis = 1)-count_per_mr_maf[t]
    
    models = {}

    for t in types:
        binom = sm.families.Binomial()
        
        lrTest = smf.glm(f'{t} + non_{t} ~ MR', 
                        count_per_mr_maf.loc[count_per_mr_maf[t]+count_per_mr_maf[f'non_{t}']>0], 
                        family=binom).fit()
        models[t]=lrTest
        fraction_per_mr_maf[f'binom_predicted_fraction_{t}'] = lrTest.predict(count_per_mr_maf[['MR']])
        models[t].save(model_out/"{t}.pickle")
    
    f, ax = plt.subplots(figsize = (3,3))
    for t, color in zip(types[:1], ['tomato', 'gold', 'seagreen', 'royalblue']):
        
        fraction_per_mr_maf.plot.scatter(x = 'MR',
                                        y = t,
                                        s = 'log_counts',
                                        ax = ax,
                                        color = color
                                        )
        
        # sub.plot(x = 'MR',
        #                 y = 'predicted_fraction_singleton', ax = ax)
        fraction_per_mr_maf.plot(x = 'MR',
                        y = f'binom_predicted_fraction_{t}', ax = ax,
                                color = color,
                                label = f'{t} model')

    plt.ylabel('Fraction')
    plt.xlabel('Mutation Rate')
    plt.legend(bbox_to_anchor=(1,0.5))
    plt.savefig(model_out/ 'fit.pdf')