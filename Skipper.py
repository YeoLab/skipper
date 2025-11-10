 ############################## SETUP #################################

# Import packages. 
import pandas as pd
from functools import reduce
import re
import os
import sys
import glob
from time import sleep
from pathlib import Path
import warnings 

# Load config file. 
locals().update(config)
WORKDIR = config.get("WORKDIR")
workdir: config['WORKDIR']
TMPDIR = config.get("TMPDIR")

# Generate directories to hold log files. 
if not os.path.exists("logs"): os.makedirs("logs")

# Set the temporary directory within the working directory by default. 
if not TMPDIR:
    config['TMPDIR'] = os.path.join(WORKDIR, "tmp")

# Check for proper overdispersion mode. 
if OVERDISPERSION_MODE not in ["clip","input"]:
    raise Exception("Overdispersion must be calculated using 'clip' or 'input' samples")

# Read and cleanup manifest.
manifest = pd.read_csv(MANIFEST, comment = "#", index_col = False).dropna(subset=['Experiment','Sample'])
manifest["CLIP_replicate"] = pd.to_numeric(manifest.CLIP_replicate, downcast="integer")
manifest["Input_replicate"] = pd.to_numeric(manifest.Input_replicate, downcast="integer")

# Remove whitespace. 
for col in manifest.columns[manifest.columns.str.contains('_fastq') | manifest.columns.str.contains('_adapter')]:
    manifest[col] = manifest[col].str.strip()

# Check that input values are valid.  
try:
    if min(manifest.groupby("Experiment")["CLIP_fastq"].agg(lambda x: len(set(x)))) < 2:
        sys.stderr.write("WARNING: NONZERO EXPERIMENTS HAVE ONLY ONE CLIP REPLICATE.\nPIPELINE MUST HALT AFTER GENERATING RAW COUNTS\nThis usually means your manifest is incorrectly formatted\n")
        print(manifest.groupby("Experiment")["CLIP_fastq"].agg(lambda x: len(set(x))))
        sleep(5)
except:
    if min(manifest.groupby("Experiment")["CLIP_fastq_1"].agg(lambda x: len(set(x)))) < 2:
        sys.stderr.write("WARNING: NONZERO EXPERIMENTS HAVE ONLY ONE CLIP REPLICATE.\nPIPELINE MUST HALT AFTER GENERATING RAW COUNTS\nThis usually means your manifest is incorrectly formatted\n")
        print(manifest.groupby("Experiment")["CLIP_fastq_1"].agg(lambda x: len(set(x))))
        sleep(5)

if max(manifest.groupby("Sample")["Input_replicate"].agg(lambda x: min(x))) > 1:
    raise Exception("Input replicates for samples in manifest do not increment from 1 as expected")

if max(manifest.groupby("Sample")["CLIP_replicate"].agg(lambda x: min(x))) > 1:
    raise Exception("CLIP replicates for samples in manifest do not increment from 1 as expected")

# create label for IN and CLIP: 
# Sample = DEK_HepG2_4020, replicate_label: DEK_HepG2_4020_IN_1 and DEK_HepG2_4020_IP_1
manifest["Input_replicate_label"] = [(str(sample) + "_IN_" + str(replicate)).replace(" ","")  for replicate, sample in zip(manifest.Input_replicate.tolist(),manifest.Sample.tolist())]
manifest["CLIP_replicate_label"] = [(str(sample) + "_IP_" + str(replicate)).replace(" ","") for replicate, sample in zip(manifest.CLIP_replicate.tolist(),manifest.Sample.tolist())]

# Extract only the relevant columns for Input and CLIP replicates (labels, fastq files, bam files, adapters).
# Drop duplicate rows to ensure each replicate is uniquely defined.
input_replicates = manifest.loc[:,manifest.columns.isin(["Input_replicate_label","Input_fastq","Input_fastq_1", "Input_fastq_2","Input_bam","Input_adapter","Input_adapter_1","Input_adapter_2"])].drop_duplicates()
clip_replicates = manifest.loc[:,manifest.columns.isin(["CLIP_replicate_label","CLIP_fastq","CLIP_fastq_1","CLIP_fastq_2","CLIP_bam","CLIP_adapter","CLIP_adapter_1","CLIP_adapter_2"])].drop_duplicates()

# Check consistency: ensure that each replicate label corresponds to exactly one set of files.
if len(input_replicates) != len(input_replicates[["Input_replicate_label"]].drop_duplicates()) or \
    len(clip_replicates) != len(clip_replicates[["CLIP_replicate_label"]].drop_duplicates()):
    raise Exception("Manifest files are not consistent across replicates")

# Collect replicate labels into lists, combine them, and store in config for downstream use.
input_replicate_labels = input_replicates.Input_replicate_label.tolist()
clip_replicate_labels = clip_replicates.CLIP_replicate_label.tolist()
replicate_labels = pd.Series(input_replicate_labels + clip_replicate_labels)
config['replicate_labels']= replicate_labels

# Map each replicate label to its corresponding FASTQ file(s) and adapter sequence(s), depending on the protocol version. 
if "Input_fastq" in manifest.columns and config['protocol']=='ENCODE4':
    config['replicate_label_to_fastqs'] = dict(zip(input_replicate_labels + clip_replicate_labels,
                                                   input_replicates.Input_fastq.tolist() + clip_replicates.CLIP_fastq.tolist()))
    config['replicate_label_to_adapter'] = dict(zip(input_replicate_labels + clip_replicate_labels,
                                                    input_replicates.Input_adapter.tolist() + clip_replicates.CLIP_adapter.tolist()))
elif config['protocol']=='ENCODE3':
    config['replicate_label_to_fastq_1'] = dict(zip(input_replicate_labels + clip_replicate_labels,
                                                    input_replicates.Input_fastq_1.tolist() + clip_replicates.CLIP_fastq_1.tolist()))
    config['replicate_label_to_fastq_2'] = dict(zip(input_replicate_labels + clip_replicate_labels,
                                                    input_replicates.Input_fastq_2.tolist() + clip_replicates.CLIP_fastq_2.tolist()))
    config['replicate_label_to_adapter_1'] = dict(zip(input_replicate_labels + clip_replicate_labels,
                                                      input_replicates.Input_adapter_1.tolist() + clip_replicates.CLIP_adapter_1.tolist()))
    config['replicate_label_to_adapter_2'] = dict(zip(input_replicate_labels + clip_replicate_labels,
                                                      input_replicates.Input_adapter_2.tolist() + clip_replicates.CLIP_adapter_2.tolist()))
else:
    raise Exception("protocol does not fit in ENCODE3 or ENCODE4")

# Map each replicate label to the expected deduplicated BAM file path. 
if config['protocol']=='ENCODE4':
    config['replicate_label_to_bams'] = dict(zip(input_replicate_labels + clip_replicate_labels, ["output/bams/dedup/genome/" + replicate_label + ".genome.Aligned.sort.dedup.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))
elif config['protocol']=='ENCODE3':
    config['replicate_label_to_bams'] = dict(zip(input_replicate_labels + clip_replicate_labels, [f"output/bams/dedup/genome_R{INFORMATIVE_READ}/" + replicate_label + f".genome.Aligned.sort.dedup.R{INFORMATIVE_READ}.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))
else:
    raise Exception("protocol does not fit in ENCODE3 or ENCODE4")

# Extract out experiment label information. 
config['experiment_labels'] = pd.Series(manifest.Experiment.drop_duplicates().tolist())
experiment_data = manifest.groupby("Experiment").agg({"CLIP_replicate_label": list, "Input_replicate_label" : list})

# Build dictionaries that link replicates together for modeling and analysis:
# Fpr determining which replicates to use when estimating variance.
config['overdispersion_replicate_lookup'] = dict(zip(manifest.CLIP_replicate_label.tolist(),
                                                     manifest.Input_replicate_label.tolist() if OVERDISPERSION_MODE == "input"
                                                     else manifest.CLIP_replicate_label.tolist()))

# For mapping each CLIP replicate label to its corresponding Input replicate label.
config['clip_to_input_replicate_label'] = dict(zip(manifest.CLIP_replicate_label.tolist(),
                                                   manifest.Input_replicate_label.tolist()))

# For each experiment, collect the set of unique Input replicate labels  and then append the CLIP replicates.
config['experiment_to_replicate_labels'] = dict(zip(experiment_data.index.tolist(),
                                                    [reduce(lambda agg, x: agg if x in agg else agg + [x], inputs, []) + clips for inputs,
                                                     clips in zip(experiment_data.Input_replicate_label, experiment_data.CLIP_replicate_label)]))

# For each experiment, map directly to its CLIP replicate labels only.
config['experiment_to_clip_replicate_labels'] = dict(zip(experiment_data.index.tolist(), experiment_data.CLIP_replicate_label))

# A nested dictionary mapping each experiment to -> input replicate to -> the *other* input replicates from the same experiment.
experiment_to_input_replicate_labels = {}
for experiment_label, label_list in zip(experiment_data.index, experiment_data.Input_replicate_label):
    # Initialize dictionary for this experiment
    experiment_to_input_replicate_labels[experiment_label] = {}
    for entry in label_list:
        replicates = set()
        # Collect all other entries except the current one
        for other_entry in label_list:
            if other_entry != entry:
                replicates.add(other_entry)
        # Map the current entry → list of its partner replicates
        experiment_to_input_replicate_labels[experiment_label].update({entry : list(replicates)})

# Save mapping into config for downstream steps
config['experiment_to_input_replicate_labels']=experiment_to_input_replicate_labels

# Add the manifest to the config file
config['manifest'] = manifest

# Collect benchmark-related outputs based on optional external mapping files.
benchmark_outputs = []

# For RBNS
if 'RBNS_MAPPING' in config:
    config['RBNS_mapping_df'] = pd.read_csv(config['RBNS_MAPPING'])
    print(config['RBNS_mapping_df'])
    experiments_to_banchmark = set(config['manifest']['Experiment']).intersection(set(config['RBNS_mapping_df']['Experiment']))
    benchmark_outputs+=[f"output/ml/benchmark/homer/RBNS/{experiment_label}.pearson_auprc.csv"
                        for experiment_label in list(experiments_to_banchmark)]
else:
    pass

# FOR SELEX
if 'SELEX_MAPPING' in config:
    config['SELEX_mapping_df'] = pd.read_csv(config['SELEX_MAPPING'])
    experiments_to_banchmark = set(config['manifest']['Experiment']).intersection(set(config['SELEX_mapping_df']['Experiment']))
    benchmark_outputs+=[f"output/ml/benchmark/homer/SELEX/{experiment_label}.pearson_auprc.csv"
                        for experiment_label in list(experiments_to_banchmark)]
else:
    pass

# Record the path of the config file that was passed to the command line (allows workflow to keep track of config file).
if '--configfile' in sys.argv:
    i = sys.argv.index('--configfile')
elif '--configfiles' in sys.argv:
    i = sys.argv.index('--configfiles')
config['CONFIG_PATH']=sys.argv[i+1]

# Make all config entries available as local variables for convenience.
locals().update(config)

# Create helper function for defining outputs of the call_enriched window rule
def call_enriched_window_output(wildcards):
    outputs = []
    for experiment_label in manifest.Experiment:
        for clip_replicate_label in config['experiment_to_clip_replicate_labels'][experiment_label]:
            outputs.append(f"output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_feature_summary.tsv")
    return outputs

# Add path to the chrom sizes file to config
config["CHROM_SIZES"] = config["STAR_DIR"] + "/chrNameLength.txt"
config["UNINFORMATIVE_READ"] = str(3 - config["INFORMATIVE_READ"])

# Used when the experiment wildcard needs to be changed to the sample label.  
experiment_to_sample = dict(zip(manifest["Experiment"], manifest["Sample"]))

############################## Define which parts of skipper to run #################################

# Always include the basic.
all_inputs = ["basic_done.txt"]

# Look through all ml related keys. 
ml_keys = ["HEADER", "RENAME", "ROULETTE_DIR", "GNOMAD_DIR", "CLINVAR_VCF", "VEP_CACHEDIR",
           "VEP_CACHE_VERSION", "SINGLETON_REFERENCE", "OE_RATIO_REFERENCE", "GNOMAD_CONSTRAINT"]

repeat_keys = ["REPEAT_TABLE", "REPEAT_BED"]

geneset_keys = ["GENE_SETS", "GENE_SET_REFERENCE", "GENE_SET_DISTANCE"]

homer_keys = ["HOMER"]

meta_keys = ["META_ANALYSIS"]

# Auto assign keys to None (avoids DAG error in snakemake).
for k in ml_keys + repeat_keys + geneset_keys + homer_keys + meta_keys:
    if k not in config or k == False:
        config[k] = ""

# A function that checks if all ml keys are present. 
def has_all_required(cfg, keys):
    return all(cfg.get(k) not in [None, ""] for k in keys)

# If all ml configs are provided, include the ml rules.
if has_all_required(config, ml_keys):
    all_inputs += [
        "ml_variants_done.txt",
        "ml_benchmark_done.txt",
        "mcross_done.txt",
    ]

if has_all_required(config, repeat_keys):
    all_inputs += [
        "repeat_done.txt",
    ]

if has_all_required(config, geneset_keys):
    all_inputs += [
        "geneset_done.txt",
    ]

if config["HOMER"] or has_all_required(config, ml_keys):
    all_inputs += [
        "homer_done.txt",
    ]

if config["META_ANALYSIS"] != "":
    all_inputs += [
        "meta_done.txt",
    ]

if config["META_ANALYSIS"] != "" and has_all_required(config, repeat_keys):
    all_inputs += [
        "meta_repeats_done.txt",
    ]

# Run rule all. 
rule all:
    input:
        all_inputs
    
############################## Call all basic outputs #################################
rule all_basic_output:
    input:
        expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam", replicate_label = replicate_labels), 
        expand("output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam.bai", replicate_label = replicate_labels), 
        expand("output/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw", replicate_label = replicate_labels),
        expand("output/bigwigs/scaled/plus/{replicate_label}.scaled.plus.bw", replicate_label = replicate_labels),
        expand("output/bigwigs/scaled/plus/{replicate_label}.scaled.cov.plus.bw", replicate_label = replicate_labels),
        expand("output/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz",
               zip, experiment_label = manifest.Experiment, clip_replicate_label = manifest.CLIP_replicate_label),
        expand("output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf", experiment_label = manifest.Experiment),
        expand("output/enrichment_reproducibility/{experiment_label}.odds_data.tsv", experiment_label = manifest.Experiment),
        lambda wildcards: call_enriched_window_output(wildcards),
        "output/figures/tsne/skipper.tsne_query.pdf",
        # Quality control
        expand("output/multiqc/{experiment_label}/multiqc_data", experiment_label = manifest.Experiment),
        expand("output/multiqc/{experiment_label}/multiqc_plots", experiment_label = manifest.Experiment),
        expand("output/multiqc/{experiment_label}/multiqc_report.html", experiment_label = manifest.Experiment),
        "output/QC/unique_fragments.csv",
        expand("output/QC/{experiment_label}.gc_bias.txt", experiment_label = manifest.Experiment),
    output:
        "basic_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

############################## Call specific optional outputs #################################

rule all_meta_output:
    input:
        expand("output/counts/genome/megatables/{genome_type}.tsv.gz", genome_type = ["feature_type_top","transcript_type_top"]),
    output:
        "meta_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

rule all_meta_repeat_output:
    input:
        expand("output/counts/repeats/megatables/{repeat_type}.tsv.gz", repeat_type = ['name', 'class', 'family']),
    output:
        "meta_repeats_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """


rule all_homer_output:
    input:
        expand("output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz", experiment_label = manifest.Experiment),
        expand("output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv", experiment_label = manifest.Experiment),
        expand("output/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.bed",experiment_label = manifest.Experiment),
        expand("output/homer/finemapped_results/{experiment_label}/homerResults.html", experiment_label = manifest.Experiment),
        expand("output/QC/{experiment_label}.nread_in_finemapped_regions.txt", experiment_label=manifest.Experiment),
    output:
        "homer_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

rule all_repeat_output:
    input:
        expand("output/counts/repeats/vectors/{replicate_label}.counts", replicate_label = replicate_labels),
        expand("output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/counts/repeats/tables/family/{experiment_label}.tsv.gz", experiment_label = manifest.Experiment),
        "output/figures/tsne_re/skipper.tsne_re_query.pdf",
    output:
        "repeat_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

rule all_geneset_output:
    input:
        expand("output/gene_sets/{experiment_label}.enriched_terms.tsv.gz", experiment_label = manifest.Experiment),
    output:
        "geneset_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

# Calls upon ouputs needed for machine learning. 
rule all_ml_variants_output:
    input:
        expand("output/ml/rbpnet_model/{experiment_label}/valid/test_data_metric.csv",
               experiment_label = manifest.Experiment),
        expand("output/ml/rbpnet_model/{experiment_label}/motif_done",
               experiment_label = manifest.Experiment),
        expand("output/variants/gnomAD_roulette/{experiment_label}.total.csv",
               experiment_label = manifest.Experiment),
        expand("output/variants/clinvar/{experiment_label}.vep.tsv",
            experiment_label = manifest.Experiment),
        expand("output/variant_analysis/{experiment_label}.clinvar_variants.csv",
               experiment_label = manifest.Experiment),
    output:
        "ml_variants_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

# Calls upon outputs needed for benchmarking. 
rule all_ml_benchmark_outputs:
    input:
        benchmark_outputs,
        expand("output/ml/benchmark/homer/{data_types}_mcross/{experiment_label}.pearson_auprc.csv",
                data_types = ['CITS'],
               experiment_label = [i for i in manifest.Experiment.tolist() if 'QKI' in i or 'RBFOX' in i or 'PUM' in i]),
        expand("output/ml/rbpnet_model_original/{experiment_label}/valid/test_data_metric.csv",
               experiment_label = [i for i in manifest.Experiment.tolist() if 'QKI' in i or 'RBFOX' in i or 'PUM' in i]),
        expand("output/ml/nt_lora/{experiment_label}/{model_name}/d_log_odds_corr.csv",
               experiment_label = [i for i in manifest.Experiment.tolist()],
                model_name = ['nucleotide-transformer-500m-human-ref']),
    output:
        "ml_benchmark_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

# Calls upon outputs needed for mcross. 
rule all_ctk:
    input:
        expand("output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.homer", experiment_label = manifest.Experiment),
        expand("output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}/{data_types}.{experiment_label}.homer",experiment_label = manifest.Experiment, data_types=['CITS']),
    output:
        "mcross_done.txt"
    resources:
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

############################## Define modules #################################

# Switch to basic
module se_preprocess:
    snakefile:
        "rules/basic/se_preprocess.smk"
    config:
        config

# Switch to basic
module pe_preprocess:
    snakefile:
        "rules/basic/pe_preprocess.smk"
    config:
        config

# Switch to basic. 
module qc:
    snakefile:
        "rules/basic/qc.smk"
    config: config

# switch to basic
module genome:
    snakefile:
        "rules/basic/genome_windows.smk"
    config: config

# switch to basic
module repeat:
    snakefile:
        "rules/basic/repeat.smk"
    config: config

# switch to basic
module finemap:
    snakefile:
        "rules/basic/finemap.smk"
    config: config

# switch to basic
module analysis:
    snakefile:
        "rules/basic/analysis.smk"
    config: config

# switch to basic
module meta_analysis:
    snakefile:
        "rules/basic/meta_analysis.smk"
    config:
        config

# Switch to basic
module bigwig:
    snakefile:
        "rules/basic/bigwig.smk"
    config:
        config

# Switch to ml
module prep_ml:
    snakefile:
        "rules/ml/prep_ml.smk"
    config:
        config

# Switch to ml
module rbpnet:
    snakefile:
        "rules/ml/train_rbpnet.smk"
    config:
        config

# switch to ml
module benchmark:
    snakefile:
        "rules/ml/benchmark_ml.smk"
    config:
        config

# switch to ml
module variants_rbpnet:
    snakefile:
        "rules/ml/variants_rbpnet.smk"
    config:
        config

# switch to mcross
module ctk_mcross:
    snakefile:
        "rules/mcross/ctk_mcross.smk"
    config:
        config

############################## Run modules #################################
    
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
use rule * from rbpnet as rbpnet_*
use rule * from variants_rbpnet as rbpnet_variants_*
use rule * from ctk_mcross
use rule * from benchmark
