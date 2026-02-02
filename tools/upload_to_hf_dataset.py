from pathlib import Path
import seqdata as sd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm
import sys

import os
from huggingface_hub import login


def sdata_to_parquet(sdata, split, exp, CHUNK=1024, outdir = '/tmp'):

    # Choose a safe chunk size
    CHUNK = 1024
    
    ds = sdata.chunk({"_sequence": CHUNK})
    
    writer = None
    
    for i in tqdm(range(0, ds.dims["_sequence"], CHUNK)):
        sub = ds.isel(_sequence=slice(i, i + CHUNK)).compute()
    
        # Convert seq bytes → strings
        seq_str = [
            b"".join(row).decode("ascii")
            for row in sub["seq"].values
        ]
    
        table = pa.Table.from_pydict({
            "chrom": sub["chrom"].astype(str).values.tolist(),
            "chromStart": sub["chromStart"].values,
            "chromEnd": sub["chromEnd"].values,
            "strand": sub["strand"].astype(str).values.tolist(),
            "name": sub["name"].values,
            "rep": sub["rep"].values,
            "gc_fraction": sub["gc_fraction"].values,
            "n_IP": sub["n_IP"].values,
            "n_IN": sub["n_IN"].values,
            "train_test": sub["train_test"].values,
            "train_val": sub["train_val"].values,
    
            # array columns
            "signal": sub["signal"].values.tolist(),
            "control": sub["control"].values.tolist(),
            "seq": seq_str,
        })
    
        if writer is None:
            writer = pq.ParquetWriter(outdir / f"{exp}.{split}.parquet", table.schema)
    
        writer.write_table(table)
    
    writer.close()

if __name__ == '__main__':
    # login
    hf_token = os.getenv("HF_TOKEN")
    login(token=hf_token, write_permission = True)

    exp = sys.argv[1]
    outdir = Path(sys.argv[2])
    data_dir = Path(f"output/ml/rbpnet_data/{exp}")
    train_sdata = sd.open_zarr(data_dir / "train.zarr").load()
    valid_sdata = sd.open_zarr(data_dir / "valid.zarr").load()
    test_sdata = sd.open_zarr(data_dir / "test.zarr").load()

    # make parquet
    sdata_to_parquet(train_sdata, exp=exp, split='train', outdir = outdir)
    sdata_to_parquet(test_sdata, exp=exp, split='test', outdir = outdir)
    sdata_to_parquet(valid_sdata, exp=exp, split='val', outdir = outdir)

    from datasets import DatasetDict, Dataset

    ds = DatasetDict({
        "train": Dataset.from_parquet(str(Path(outdir)/f"{exp}.train.parquet")),
        "validation": Dataset.from_parquet(str(Path(outdir)/f"{exp}.val.parquet")),
        "test": Dataset.from_parquet(str(Path(outdir)/f"{exp}.test.parquet")),
    })

    ds.push_to_hub(
        "yeolab/eCLIP",
        config_name=exp
    )
