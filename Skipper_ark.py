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
GINI_CUTOFF = config.get("GINI_CUTOFF")
BLACKLIST = config.get("BLACKLIST")
NORMALIZATION_MODE = config.get("NORMALIZATION_MODE")
THRESHOLD_MIN = config.get("THRESHOLD_MIN")
HOMER = config.get("HOMER")

############################ ESTABLISH DEFAULTS ###########################

# Set the temporary directory within the working directory by default. 
if not TMPDIR:
    config['TMPDIR'] = os.path.join(WORKDIR, "tmp")

if not GINI_CUTOFF:
    config['GINI_CUTOFF'] = 0.9

if not BLACKLIST or str(BLACKLIST).strip().lower() in {"none", "null", "na", "n/a"}:
    BLACKLIST = None

if not HOMER:
    HOMER = ""

if not NORMALIZATION_MODE:
    NORMALIZATION_MODE = "new"

if not THRESHOLD_MIN:
    THRESHOLD_MIN = 2

config["BLACKLIST"] = BLACKLIST
config["HOMER"] = HOMER
config["NORMALIZATION_MODE"] = NORMALIZATION_MODE
config["THRESHOLD_MIN"] = THRESHOLD_MIN

# Check for proper overdispersion mode. 
if OVERDISPERSION_MODE not in ["clip","input"]:
    raise Exception("Overdispersion must be calculated using 'clip' or 'input' samples")

########################## CLEAN MANIFEST #############################

# Read and cleanup manifest.
manifest = pd.read_csv(MANIFEST, comment = "#", index_col = False).dropna(subset=['Experiment','Sample'])
manifest["CLIP_replicate"] = pd.to_numeric(manifest.CLIP_replicate, downcast="integer")
manifest["Input_replicate"] = pd.to_numeric(manifest.Input_replicate, downcast="integer")

# Remove whitespace. 
for col in manifest.columns[manifest.columns.str.contains('_fastq') | manifest.columns.str.contains('_adapter')]:
    manifest[col] = manifest[col].str.strip()

# Check that input values are valid.
def _first_present(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

clip_id_col = _first_present(manifest, ["CLIP_bam", "CLIP_fastq", "CLIP_fastq_1"])
if clip_id_col is None:
    raise Exception("Manifest must contain one of: CLIP_bam, CLIP_fastq, CLIP_fastq_1")

clip_reps_per_experiment = manifest.groupby("Experiment")[clip_id_col].agg(lambda x: len(set(x)))
if clip_reps_per_experiment.min() < 2:
    sys.stderr.write(
        "WARNING: NONZERO EXPERIMENTS HAVE ONLY ONE CLIP REPLICATE.\n"
        "PIPELINE MUST HALT AFTER GENERATING RAW COUNTS\n"
        "This usually means your manifest is incorrectly formatted\n"
    )
    print(clip_reps_per_experiment)
    sleep(5)

# Additional index checks.
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
if "CLIP_bam" not in manifest.columns:
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
if "CLIP_bam" in manifest.columns:
    config['replicate_label_to_bams'] = dict(zip(input_replicate_labels + clip_replicate_labels, input_replicates.Input_bam.to_list() + clip_replicates.CLIP_bam.to_list()))
elif config['protocol']=='ENCODE4':
    config['replicate_label_to_bams'] = dict(zip(input_replicate_labels + clip_replicate_labels, ["output/secondary_results/bams/dedup/genome/" + replicate_label + ".genome.Aligned.sort.dedup.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))
elif config['protocol']=='ENCODE3':
    config['replicate_label_to_bams'] = dict(zip(input_replicate_labels + clip_replicate_labels, [f"output/secondary_results/bams/dedup/genome_R{INFORMATIVE_READ}/" + replicate_label + f".genome.Aligned.sort.dedup.R{INFORMATIVE_READ}.bam" for replicate_label in input_replicate_labels + clip_replicate_labels] ))
else:
    raise Exception("pre-processing protocol does not fit in ENCODE3 or ENCODE4, and no pre-processed bams provided")

# Extract out experiment label information. 
config['experiment_labels'] = pd.Series(manifest.Experiment.drop_duplicates().tolist())
experiment_data = manifest.groupby("Experiment").agg({"CLIP_replicate_label": list, "Input_replicate_label" : list})

########################## Build dictionaries that link replicates together for modeling and analysis ############################

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

########################## FINAL SETUP ###########################

# Add the manifest to the config file
config['manifest'] = manifest

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
            outputs.append(f"output/secondary_results/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_feature_summary.tsv")
    return outputs

# Add path to the chrom sizes file to config
config["CHROM_SIZES"] = config["STAR_DIR"] + "/chrNameLength.txt"
config["UNINFORMATIVE_READ"] = str(3 - config["INFORMATIVE_READ"])

############################## Define which parts of skipper to run #################################

# Always include the basic.
all_inputs = ["basic_done.txt"]

repeat_keys = ["REPEAT_TABLE", "REPEAT_BED"]

geneset_keys = ["GENE_SETS", "GENE_SET_REFERENCE", "GENE_SET_DISTANCE"]

homer_keys = ["HOMER"]

# Auto assign keys to None (avoids DAG error in snakemake).
for k in repeat_keys + geneset_keys + homer_keys:
    if k not in config or k == False:
        config[k] = "DO/NOT/RUN"

# A function that checks if all keys are present. 
def has_all_required(cfg, keys):
    return all(cfg.get(k) not in [None, "DO/NOT/RUN"] for k in keys)

if "CLIP_bam" not in manifest.columns:
    all_inputs += [
        "QC_done.txt"
    ]

if has_all_required(config, repeat_keys):
    all_inputs += [
        "repeat_done.txt",
    ]

if has_all_required(config, geneset_keys):
    all_inputs += [
        "geneset_done.txt",
    ]

if config["HOMER"] != "":
    all_inputs += [
        "homer_done.txt",
    ]

print(all_inputs)

# Run rule all. 
rule all:
    input:
        all_inputs
    
############################## Call all basic outputs #################################
rule all_basic_output:
    input:
        list(config["replicate_label_to_bams"].values()),
        [p + ".bai" for p in list(config["replicate_label_to_bams"].values())],
        expand("output/secondary_results/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw", replicate_label=replicate_labels),
        expand("output/secondary_results/bigwigs/scaled/plus/{replicate_label}.scaled.plus.bw", replicate_label = replicate_labels),
        expand("output/secondary_results/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz",
               zip, experiment_label = manifest.Experiment, clip_replicate_label = manifest.CLIP_replicate_label),
        expand("output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/figures/secondary_figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf", experiment_label = manifest.Experiment),
        expand("output/secondary_results/enrichment_reproducibility/{experiment_label}.odds_data.tsv", experiment_label = manifest.Experiment),
        lambda wildcards: call_enriched_window_output(wildcards),
        "output/figures/tsne/skipper.tsne_query.pdf",
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

rule all_qc:
    input:
        expand("output/QC/multiqc/{experiment_label}/multiqc_data", experiment_label = manifest.Experiment),
        expand("output/QC/multiqc/{experiment_label}/multiqc_plots", experiment_label = manifest.Experiment),
        expand("output/QC/multiqc/{experiment_label}/multiqc_report.html", experiment_label = manifest.Experiment),
    output:
        "QC_done.txt"
    resources: 
        mem_mb=400,
        run_time=20
    shell:
        """
        touch {output}
        """

rule all_homer_output:
    input:
        expand("output/secondary_results/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz", experiment_label = manifest.Experiment),
        expand("output/secondary_results/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv", experiment_label = manifest.Experiment),
        expand("output/secondary_results/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.bed",experiment_label = manifest.Experiment),
        expand("output/homer/finemapped_results/{experiment_label}/homerResults.html", experiment_label = manifest.Experiment),
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
        expand("output/secondary_results/counts/repeats/vectors/{replicate_label}.counts", replicate_label = replicate_labels),
        expand("output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz", experiment_label = manifest.Experiment),
        expand("output/secondary_results/counts/repeats/tables/family/{experiment_label}.tsv.gz", experiment_label = manifest.Experiment),
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

############################## Define modules #################################


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
    config:
        config


module genome:
    snakefile:
        "rules/genome_windows.smk"
    config: 
        config


module repeat:
    snakefile:
        "rules/repeat.smk"
    config: 
        config


module finemap:
    snakefile:
        "rules/finemap.smk"
    config: 
        config


module analysis:
    snakefile:
        "rules/analysis.smk"
    config:
        config

module bigwig:
    snakefile:
        "rules/bigwig.smk"
    config:
        config


############################## Run modules #################################

if "CLIP_bam" in manifest.columns:
    pass
elif config['protocol']=='ENCODE4':
    use rule * from se_preprocess as se_*
else:
    use rule * from pe_preprocess as pe_*

use rule * from bigwig
use rule * from qc
use rule * from genome
use rule * from repeat
use rule * from finemap
use rule consult_encode_reference_re from analysis
use rule consult_encode_reference from analysis
if has_all_required(config, geneset_keys):
    use rule consult_term_reference from analysis
if config["HOMER"] != "":
    use rule sample_background_windows_by_region from analysis
    use rule run_homer from analysis
