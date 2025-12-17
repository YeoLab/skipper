import pandas as pd
import sys

if __name__=='__main__':
    input = sys.argv[1]
    feature_output = sys.argv[2]
    transcript_output = sys.argv[3]


    fcounts = []
    tcounts = []
    for df in pd.read_csv(input, sep = '\t',
                        chunksize = 1000):
        sample_col = df.columns[17:]
        fcount = df.groupby(by = 'feature_type_top')[sample_col].sum().reset_index()
        tcount = df.groupby(by = 'transcript_type_top')[sample_col].sum().reset_index()
        
        
        fcounts.append(fcount)
        tcounts.append(tcount)
    fcounts = pd.concat(fcounts, axis = 0).groupby(by = 'feature_type_top').sum()
    tcounts = pd.concat(tcounts, axis = 0).groupby(by = 'transcript_type_top').sum()

    fcounts.to_csv(feature_output, sep = '\t')
    tcounts.to_csv(transcript_output, sep = '\t')