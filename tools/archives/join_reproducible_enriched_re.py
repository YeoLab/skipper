from pathlib import Path
import pandas as pd
import sys

def reproducible_re(input_path):
    ''' given Skipper output/ path, count how many reproducible enriched windows per feature_type and transcript_type '''
    families = [] # save full data
    families_l2or = []
    
    for f in (Path(input_path) / 'output/reproducible_enriched_re').glob('*tsv.gz'):
        df = pd.read_csv(f, sep = '\t', index_col = 0)
        
        # append binary matrix
        name = f.name.split('.')[0]
        df[name]=True
        families.append(df[name])

        # append l2or matrix
        l2or = df['enrichment_l2or_mean']
        l2or.name = name
        families_l2or.append(l2or)
        
    
    families = pd.concat(families, axis = 1).fillna(False)
    families_l2or = pd.concat(families_l2or, axis = 1)

    return families, families_l2or


if __name__ == '__main__':
    input_path = sys.argv[1]
    binary_output = sys.argv[2]
    l2or_output = sys.argv[3]

    binary_matrix, l2or_matrix = reproducible_re(input_path)

    binary_matrix.to_csv(binary_output)
    l2or_matrix.to_csv(l2or_output)