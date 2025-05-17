import pandas as pd
import sys
from pathlib import Path

if __name__ == '__main__':
    exp = sys.argv[1]
    indir = Path(sys.argv[2])
    outf = Path(sys.argv[3])

    def annotate_score(chr, exp):
        print('processing chr' + str(chr))
        variants = pd.read_csv(indir / f'{chr}.vcf', 
                        sep = '\t', names =['CHROM','POS','ID','REF','ALT','INFO/AC','INFO/AN','INFO/MR','INFO/AR','INFO/MG','INFO/MC'],
            na_values = '.')
        score = pd.read_csv(indir / f'{exp}.{chr}.score.csv',
                    index_col = 0)
        score[['CHROM', 'POS', 'NU', 'name']]=score['ID'].str.split('-', expand = True)
        score['POS']=score['POS'].astype(int)
        ref_score = score.merge(variants, left_on = ['CHROM', 'POS', 'NU'],
                right_on = ['CHROM', 'POS', 'REF'],
                        how = 'right')[['CHROM', 'POS', 'REF', 'ALT', 'dlogodds_pred']]
        alt_score = score.merge(variants, left_on = ['CHROM', 'POS', 'NU'],
                    right_on = ['CHROM', 'POS', 'ALT'],
                            how = 'right')[['CHROM', 'POS', 'REF', 'ALT', 'dlogodds_pred','INFO/AC', 'INFO/AN', 'INFO/MR',
            'INFO/AR', 'INFO/MG', 'INFO/MC', 'name']]
        all_scores = ref_score.merge(alt_score, 
                                    left_on = ['CHROM', 'POS', 'REF', 'ALT'],
                                    right_on = ['CHROM', 'POS','REF', 'ALT'],
                                    suffixes = ('_REF', '_ALT')
                                )
        all_scores['delta_score'] = all_scores['dlogodds_pred_ALT']-all_scores['dlogodds_pred_REF']
        return all_scores
    all_scores = pd.concat([annotate_score(c, exp) for c in range(1,23)],
                       axis = 0)
    all_scores.to_csv(outf, index = False)