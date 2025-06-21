import seqdata as sd
import pandas as pd
import os
from scipy.stats import pearsonr
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
from pathlib import Path
import sys
import torch
import numpy as np

if __name__ == '__main__':


    skipper_dir = Path(sys.argv[1])
    exp = sys.argv[2]
    motif_source = sys.argv[3]
    rbpnet_path = sys.argv[4]
    outf = sys.argv[5]
    sys.path.append(rbpnet_path)
    from metrics import dlog_odds_from_data

    data_dir = skipper_dir / "output/ml/rbpnet_data" / exp
    test_sdata = sd.open_zarr(os.path.join(data_dir, "test.zarr")).load()


    reproducible_enriched_windows = pd.read_csv(
        skipper_dir
        / f"output/reproducible_enriched_windows/{exp}.reproducible_enriched_windows.tsv.gz",
        sep="\t",
    )

    
    homer_predictions = pd.read_csv(skipper_dir / f"output/ml/benchmark/homer/{motif_source}/{exp}.csv",
                                sep = '\t',
                                names = ['name', 'position', 'sequence', 'motif', 'strand', 'score'])
    homer_predictions['name']=homer_predictions['name'].str.split('-', expand = True)[0]
    homer_predictions.drop_duplicates(inplace = True)
    homer_predictions.dropna(subset = ['name'], inplace = True)
    homer_predictions['name'] = homer_predictions['name'].astype(int)

    per_window_homer_score = homer_predictions.groupby(by = ['name', 'motif'])['score'].sum().to_frame().reset_index()
    per_window_homer_score['is_reproducible_enriched_window'] = per_window_homer_score.index.isin(reproducible_enriched_windows['name'].tolist())

    y_logodd, y_dlogodd = dlog_odds_from_data(
    {
        "n_IP": torch.from_numpy(test_sdata["n_IP"].values),
        "n_IN": torch.from_numpy(test_sdata["n_IN"].values),
        "gc_fraction": torch.from_numpy(test_sdata["gc_fraction"].values),
    }
    )

    test_metrics = pd.DataFrame(
        {
            "dlogodds": y_dlogodd,
            "logodds": y_logodd,
            "n_IP": torch.from_numpy(test_sdata["n_IP"].values),
            "n_IN": torch.from_numpy(test_sdata["n_IN"].values),
            "name": test_sdata["name"].values,
            "strand": test_sdata["strand"].values,
        }
    )
    test_metrics["total"] = test_metrics["n_IP"] + test_metrics["n_IN"]
    per_window_homer_score = per_window_homer_score.merge(test_metrics.reset_index(), on = 'name')
    data = []
    for name, group in per_window_homer_score.groupby(by = 'motif'):
        precision, recall, thresholds = precision_recall_curve(
        group['is_reproducible_enriched_window'], group['score'].astype(float)
        )
        # Use AUC function to calculate the area under the curve of precision recall curve
        auprc = auc(recall, precision)
        
        try:
            r, pval = pearsonr(group['score'], group['dlogodds'])
        except Exception as e:
            print('all', e)
            r = np.nan
            
        scores = pd.Series({'AUPRC':auprc,
                'pearsonr': r,
                            'motif': name,
                            })
                
        for nread in [10, 50, 100, 200]:
            try:
                sub = group[group["total"] > nread]
                r, p = pearsonr(sub["dlogodds"], sub["score"])
                scores[f'dlogodds pearson(total>{nread})']=r
            except Exception as e:
                print(nread, e)
        data.append(scores)
    pd.concat(data, axis = 1).T.to_csv(outf)