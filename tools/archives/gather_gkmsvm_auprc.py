import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve, auc
import matplotlib.pyplot as plt
from pathlib import Path
import sys
def mean_auprc(df):
    aucs = []
    for name, group in df.groupby(by = 'CV-set'):
        precision, recall, thresholds = precision_recall_curve(group['label'], group['SVM score'])

        aucs.append(auc(recall, precision))

    mean_auprc, std_auprc = np.mean(aucs), np.std(aucs)
    return mean_auprc, std_auprc

if __name__=='__main__':
    output = sys.argv[1]

    data = []
    for f in Path('output/ml/gkmsvm/').glob('*cvpred.txt'):
        df = pd.read_csv(f, 
                sep = '\t',
            names = ['sequence_id', 'SVM score', 'label', 'CV-set'])
        mean, std = mean_auprc(df)
        
        data.append([f.name.split('.')[0],mean, std])
    data = pd.DataFrame(data, columns = ['Experiment', 'mean AUPRC', 'std AUPRC'])

    data.to_csv(output)