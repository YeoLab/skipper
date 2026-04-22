import numpy as np
import matplotlib.pyplot as plt
import vizsequence
import sys
from pathlib import Path



if __name__ == '__main__':
    in_file = sys.argv[1]
    outdir = Path(sys.argv[2])
    impscores = [
        np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
        for x in open(in_file)
        ]



    #visualize importance scores on a couple of sequences, as a sanity check

    for i,score in enumerate(impscores):
        vizsequence.viz_sequence.plot_weights(score)
        plt.savefig(outdir/f'{i}.pdf')
        
