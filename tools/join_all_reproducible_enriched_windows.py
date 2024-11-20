from pathlib import Path
import pandas as pd
import sys
from pathlib import Path

if __name__ == '__main__':
    indir = Path(sys.argv[1])
    outf = sys.argv[2]

    all_data = []

    for f in (indir/'output/reproducible_enriched_windows').glob('*.tsv.gz'):
        rbp = f.name.split('.')[0]
        reproducible_windows = pd.read_csv(f, sep = '\t')

        tested_windows_1 = pd.read_csv(indir / 'output/tested_windows'/ f'{rbp}.{rbp}_IP_1.tested_windows.tsv.gz', sep = '\t')
        tested_windows_2 = pd.read_csv(indir / 'output/tested_windows'/ f'{rbp}.{rbp}_IP_2.tested_windows.tsv.gz', sep = '\t')

        both_tested = set(tested_windows_1['name']).intersection(set(tested_windows_2['name']))

        value = pd.Series(index = both_tested)
        value[reproducible_windows['name']]=True
        value.fillna(False, inplace = True)
        value.name = rbp
        all_data.append(value)

        
    all_data = pd.concat(all_data, axis = 1)
    all_data.to_csv(outf)