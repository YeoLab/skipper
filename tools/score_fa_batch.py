"""Score many fasta files with one RBPNet model, loading the model once.

usage: python score_fa_batch.py <model_dir> <zarr_outdir> <fa_1> <out_1> [<fa_2> <out_2> ...]
Requires RBPNet (RBPNET_PATH) on PYTHONPATH. Empty fasta files get empty outputs,
matching score_fa in rules/variants_rbpnet.smk.
"""
import sys
import os
from pathlib import Path
from RBPNet import *
from inference import *
from module import Module

if __name__ == '__main__':
    log_dir = Path(sys.argv[1])
    zarr_outdir = Path(sys.argv[2])
    pairs = sys.argv[3:]
    assert len(pairs) % 2 == 0, 'expected <fa> <out> pairs'
    zarr_outdir.mkdir(parents=True, exist_ok=True)

    module = Module.load_from_checkpoint(find_latest_ckp(log_dir),
                                         arch = RBPNet(mask = 100),
                                         )
    for fa, outf in zip(pairs[::2], pairs[1::2]):
        if os.path.getsize(fa) == 0:
            open(outf, 'w').close()
            continue
        score_fa(Path(fa), module, zarr_outdir = zarr_outdir).to_csv(outf)
