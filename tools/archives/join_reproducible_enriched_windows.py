from pathlib import Path
import pandas as pd
import sys
def reproducible_windows(input_path, feature_type):
    ''' given Skipper output/ path, count how many reproducible enriched windows per feature_type and transcript_type '''
    families = [] # save full data
    families_l2or = []
    
    for f in (Path(input_path) / 'output/reproducible_enriched_windows').glob('*tsv.gz'):
        df = pd.read_csv(f, sep = '\t')
        df = df.loc[df['feature_type_top'].isin(feature_type)]
        
        # join binary
        name = f.name.split('.')[0]
        df[name]=True
        families.append(df.set_index('name')[name])

        # join l2or
        try:
            l2or = df.set_index('name')['enrichment_l2or_mean']
            l2or.name = name
            families_l2or.append(l2or)
        except KeyError as e:
            # some don't have any reproducible enrichment
            pass

        
    
    families = pd.concat(families, axis = 1).fillna(False)
    families_l2or = pd.concat(families_l2or, axis = 1)
    return families, families_l2or

if __name__=='__main__':
    input_path = sys.argv[1]
    feature_type = sys.argv[2].split(',')
    binary_output = sys.argv[3]
    l2or_output = sys.argv[4]

    binary_matrix, l2or_matrix = reproducible_windows(input_path, feature_type)

    binary_matrix.to_csv(binary_output)
    l2or_matrix.to_csv(l2or_output)
                           