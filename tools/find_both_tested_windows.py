import sys
import pandas as pd
from pybedtools import BedTool

if __name__ == '__main__':
    tested_windows = pd.concat([pd.read_csv(f, sep = '\t') for f in sys.argv[1].split(' ')],
                               axis = 0)
    tested_times = tested_windows['name'].value_counts()
    both_tested = tested_times[tested_times >= 2].index

    cols = ['chr', 'start', 'end', 'name', 'gc', 'strand']
    annot = tested_windows[cols].drop_duplicates()
    both_tested_windows = annot.loc[annot['name'].isin(both_tested)]
    both_tested_windows.to_csv(sys.argv[2], sep = '\t', index = False, header = False)
    print(both_tested_windows.shape)

    bed = BedTool.from_dataframe(both_tested_windows)
    # merge adjacent ones
    bed.sort().merge(c = [4,5,6], s = True, o = ['distinct', 'mean', 'distinct']).saveas(sys.argv[3])
