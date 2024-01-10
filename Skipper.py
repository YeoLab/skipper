import pandas as pd
from functools import reduce
import re
import os
import sys
import glob
from time import sleep
from pathlib import Path
import warnings 

# example command
"""
# ENCODE 4
snakemake -kps Skipper.py \
    -j 30 \
    --cluster "qsub -e {params.error_file} -o {params.out_file} -l walltime={params.run_time} -l nodes=1:ppn={threads} -q home-yeo" \
    --configfile /home/hsher/projects/skipper/encore_configs/Skipper_config_small_test.yaml \
    --conda-prefix /home/hsher/snakeconda \
    --use-conda \
    --use-singularity \
    --singularity-prefix /home/hsher/scratch/singularity \
    --singularity-args "--bind /oasis --bind /projects --bind /scratch" \
    -n

# ENCODE 3
snakemake -kps Skipper.py \
    -j 30 \
    --cluster "qsub -e {params.error_file} -o {params.out_file} -l walltime={params.run_time} -l nodes=1:ppn={threads} -q home-yeo" \
    --configfile /home/hsher/projects/skipper/encode_configs/Skipper_pe_small_test.yaml \
    --conda-prefix /home/hsher/snakeconda \
    --use-conda \
    --use-singularity \
    --singularity-prefix /home/hsher/scratch/singularity \
    --singularity-args "--bind /oasis --bind /projects --bind /scratch" \
    output/ml/gkmsvm/RBFOX2_HepG2_ENCSR987FTF.cvpred.txt


"""

locals().update(config)

workdir: config['WORKDIR']

if not os.path.exists("stderr"): os.makedirs("stderr")
if not os.path.exists("stdout"): os.makedirs("stdout")


if OVERDISPERSION_MODE not in ["clip","input"]:
    raise Exception("Overdispersion must be calculated using 'clip' or 'input' samples")

# read and cleanup manifest
manifest = pd.read_csv(MANIFEST, comment = "#", index_col = False).dropna(subset=['Experiment','Sample'])
manifest["CLIP_replicate"] = pd.to_numeric(manifest.CLIP_replicate, downcast="integer")
manifest["Input_replicate"] = pd.to_numeric(manifest.Input_replicate, downcast="integer")

for col in manifest.columns[manifest.columns.str.contains('_fastq') | manifest.columns.str.contains('_adapter')]:
    manifest[col] = manifest[col].str.strip()

# check numbers are correct
try:
    if min(manifest.groupby("Experiment")["CLIP_fastq"].agg(lambda x: len(set(x)))) < 2:
        sys.stderr.write("WARNING: NONZERO EXPERIMENTS HAVE ONLY ONE CLIP REPLICATE.\nPIPELINE MUST HALT AFTER GENERATING RAW COUNTS\nThis usually means your manifest is incorrectly formatted\n")
        sleep(5)
except:
    if min(manifest.groupby("Experiment")["CLIP_fastq_1"].agg(lambda x: len(set(x)))) < 2:
        sys.stderr.write("WARNING: NONZERO EXPERIMENTS HAVE ONLY ONE CLIP REPLICATE.\nPIPELINE MUST HALT AFTER GENERATING RAW COUNTS\nThis usually means your manifest is incorrectly formatted\n")
        sleep(5)

if max(manifest.groupby("Sample")["Input_replicate"].agg(lambda x: min(x))) > 1:
    raise Exception("Input replicates for samples in manifest do not increment from 1 as expected")

if max(manifest.groupby("Sample")["CLIP_replicate"].agg(lambda x: min(x))) > 1:
    raise Exception("CLIP replicates for samples in manifest do not increment from 1 as expected")

# create label for IN and CLIP: 
# Sample = DEK_HepG2_4020, replicate_label: DEK_HepG2_4020_IN_1 and DEK_HepG2_4020_IP_1
manifest["Input_replicate_label"] = [(str(sample) + "_IN_" + str(replicate)).replace(" ","")  for replicate, sample in zip(manifest.Input_replicate.tolist(),manifest.Sample.tolist())]
manifest["CLIP_replicate_label"] = [(str(sample) + "_IP_" + str(replicate)).replace(" ","") for replicate, sample in zip(manifest.CLIP_replicate.tolist(),manifest.Sample.tolist())]

input_replicates = manifest.loc[:,manifest.columns.isin(["Input_replicate_label","Input_fastq","Input_fastq_1", "Input_fastq_2","Input_bam","Input_adapter","Input_adapter_1","Input_adapter_2"])].drop_duplicates()
clip_replicates = manifest.loc[:,manifest.columns.isin(["CLIP_replicate_label","CLIP_fastq","CLIP_fastq_1","CLIP_fastq_2","CLIP_bam","CLIP_adapter","CLIP_adapter_1","CLIP_adapter_2"])].drop_duplicates()



if len(input_replicates) != len(input_replicates[["Input_replicate_label"]].drop_duplicates()) or \
    len(clip_replicates) != len(clip_replicates[["CLIP_replicate_label"]].drop_duplicates()):
    raise Exception("Manifest files are not consistent across replicates")

input_replicate_labels = input_replicates.Input_replicate_label.tolist()
clip_replicate_labels = clip_replicates.CLIP_replicate_label.tolist()
replicate_labels = pd.Series(input_replicate_labels + clip_replicate_labels)
config['replicate_labels']= replicate_labels

# if all(bam in manifest.columns.tolist() for bam in ["Input_bam", "CLIP_bam"]):
#     replicate_label_to_bams = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_bam.tolist() + clip_replicates.CLIP_bam.tolist()))    
# else:

# FASTQ and ADAPTOR
if "Input_fastq" in manifest.columns and config['protocol']=='ENCODE4':
    config['replicate_label_to_fastqs'] = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_fastq.tolist() + clip_replicates.CLIP_fastq.tolist()))
    config['replicate_label_to_adapter'] = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_adapter.tolist() + clip_replicates.CLIP_adapter.tolist()))
elif config['protocol']=='ENCODE3':
    config['replicate_label_to_fastq_1'] = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_fastq_1.tolist() + clip_replicates.CLIP_fastq_1.tolist()))
    config['replicate_label_to_fastq_2'] = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_fastq_2.tolist() + clip_replicates.CLIP_fastq_2.tolist()))
    config['replicate_label_to_adapter_1'] = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_adapter_1.tolist() + clip_replicates.CLIP_adapter_1.tolist()))
    config['replicate_label_to_adapter_2'] = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_adapter_2.tolist() + clip_replicates.CLIP_adapter_2.tolist()))
else:
    raise Exception("protocol does not fit in ENCODE3 or ENCODE4")

# BAMs
if config['protocol']=='ENCODE4':
    config['replicate_label_to_bams'] = dict(zip(input_replicate_labels + clip_replicate_labels, ["output/bams/dedup/genome/" + replicate_label + ".genome.Aligned.sort.dedup.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))
elif config['protocol']=='ENCODE3':
    config['replicate_label_to_bams'] = dict(zip(input_replicate_labels + clip_replicate_labels, [f"output/bams/dedup/genome_R{INFORMATIVE_READ}/" + replicate_label + f".genome.Aligned.sort.dedup.R{INFORMATIVE_READ}.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))
else:
    raise Exception("protocol does not fit in ENCODE3 or ENCODE4")

# EXPERIMENT LABELS
config['experiment_labels'] = pd.Series(manifest.Experiment.drop_duplicates().tolist())
experiment_data = manifest.groupby("Experiment").agg({"CLIP_replicate_label": list, "Input_replicate_label" : list})


# OVERDISPERSION and BACKGROUND PAIRING
config['overdispersion_replicate_lookup'] = dict(zip(manifest.CLIP_replicate_label.tolist(), manifest.Input_replicate_label.tolist() if OVERDISPERSION_MODE == "input" else manifest.CLIP_replicate_label.tolist()))
config['clip_to_input_replicate_label'] = dict(zip(manifest.CLIP_replicate_label.tolist(), manifest.Input_replicate_label.tolist()))
config['experiment_to_replicate_labels'] = dict(zip(experiment_data.index.tolist(), [reduce(lambda agg, x: agg if x in agg else agg + [x], inputs, []) + clips for inputs, clips in zip(experiment_data.Input_replicate_label, experiment_data.CLIP_replicate_label)]))
config['experiment_to_clip_replicate_labels'] = dict(zip(experiment_data.index.tolist(), experiment_data.CLIP_replicate_label))

experiment_to_input_replicate_labels = {}
for experiment_label, label_list in zip(experiment_data.index, experiment_data.Input_replicate_label):
    experiment_to_input_replicate_labels[experiment_label] = {}
    for entry in label_list:
        replicates = set()
        for other_entry in label_list:
            if other_entry != entry:
                replicates.add(other_entry)
        experiment_to_input_replicate_labels[experiment_label].update({entry : list(replicates)})
config['experiment_to_input_replicate_labels']=experiment_to_input_replicate_labels

# Fool-proof Detect misalignment for GFF and PARTITION
if Path(GFF).name.replace('.gff3.gz', '') != Path(FEATURE_ANNOTATIONS).name.replace('.tiled_partition.features.tsv.gz', ''):
    warnings.warn(f'''Detected Name Mismatch in GFF and FEATURE ANNOTATIONS:
    FEATURE_ANNOTATIONS={FEATURE_ANNOTATIONS}
    GFF={GFF}
    Check if they are the same cell line
    ''')

config['manifest'] = manifest


rule all:
    input:
        # expand("output/fastqc/initial/{replicate_label}_fastqc.html", replicate_label = replicate_labels
        # ) if config['protocol']=='ENCODE4' else ,
        # expand("output/fastqc/processed/{replicate_label}.trimmed.umi_fastqc.html", replicate_label = replicate_labels), 
        expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam", replicate_label = replicate_labels), 
        expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam.bai", replicate_label = replicate_labels), 
        expand("output/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw", replicate_label = replicate_labels),
        expand("output/bigwigs/scaled/plus/{replicate_label}.scaled.plus.bw", replicate_label = replicate_labels),
        expand("output/counts/repeats/vectors/{replicate_label}.counts", replicate_label = replicate_labels),
        expand("output/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz", zip, experiment_label = manifest.Experiment, clip_replicate_label = manifest.CLIP_replicate_label),
        expand("output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf", experiment_label = manifest.Experiment),
        expand("output/counts/repeats/tables/family/{experiment_label}.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz", experiment_label = manifest.Experiment),
        expand("output/homer/finemapped_results/{experiment_label}/homerResults.html", experiment_label = manifest.Experiment),
        expand("output/gene_sets/{experiment_label}.enriched_terms.tsv.gz", experiment_label = manifest.Experiment),
        "output/figures/tsne/skipper.tsne_query.pdf",
        expand("output/multiqc/{experiment_label}/multiqc_data", experiment_label = manifest.Experiment),
        expand("output/multiqc/{experiment_label}/multiqc_plots", experiment_label = manifest.Experiment),
        expand("output/multiqc/{experiment_label}/multiqc_report.html", experiment_label = manifest.Experiment),
        expand("output/counts/genome/megatables/{genome_type}.tsv.gz", genome_type = ["feature_type_top","transcript_type_top"]),
        expand("output/counts/repeats/megatables/{repeat_type}.tsv.gz", repeat_type = ['name', 'class', 'family']),
        "output/QC/unique_fragments.csv",
        # expand("output/ml/sequence/{experiment_label}.foreground.fa", experiment_label = manifest.Experiment),
        # expand("output/ml/gkmsvm/{experiment_label}.cvpred.txt", experiment_label = manifest.Experiment),
        # expand("output/ml/gkmsvm/{experiment_label}.model.txt", experiment_label = manifest.Experiment),
        # "output/ml/gkmsvm/AUPRC.txt",
        # "variants_done.txt"
module se_preprocess:
    snakefile:
        "rules/se_preprocess.smk"
    config:
        config

module pe_preprocess:
    snakefile:
        "rules/pe_preprocess.smk"
    config:
        config

module qc:
    snakefile:
        "rules/qc.smk"
    config: config

module genome:
    snakefile:
        "rules/genome_windows.smk"
    config: config

module repeat:
    snakefile:
        "rules/repeat.smk"
    config: config

module finemap:
    snakefile:
        "rules/finemap.smk"
    config: config

module analysis:
    snakefile:
        "rules/analysis.smk"
    config: config

module meta_analysis:
    snakefile:
        "rules/meta_analysis.smk"
    config:
        config
module bigwig:
    snakefile:
        "rules/bigwig.smk"
    config:
        config
module prep_ml:
    snakefile:
        "rules/prep_ml.smk"
    config:
        config

module variants:
    snakefile:
        "rules/variants.smk"
    config:
        config



if config['protocol']=='ENCODE4':
    use rule * from se_preprocess as se_*
else:
    use rule * from pe_preprocess as pe_*

use rule * from bigwig
use rule * from qc
use rule * from genome
use rule * from repeat
use rule * from finemap
use rule * from analysis
use rule * from meta_analysis
use rule * from prep_ml as ml_*
use rule * from variants as variants_*

