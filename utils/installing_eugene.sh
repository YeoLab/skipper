conda create -n eugene2 

conda activate eugene2
mamba install python==3.9.2
python -m pip install eugene-tools pandas seaborn seqdata statsmodels seqpro torch scipy numpy tqdm pybedtools xarray tensorboard tabulate
python -m pip install eugene-tools==0.1.2 pandas seaborn seqdata==0.1.3 statsmodels seqpro==0.1.11 torch scipy numpy tqdm pybedtools xarray tensorboard tabulate modisco-lite==2.2.1