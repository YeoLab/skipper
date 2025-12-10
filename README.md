# skipper
![Skipper cartoon](documents/logo.png)

Skip the peaks and expose RNA-binding in CLIP data

See published article in Cell Genomics: https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00085-X

## Yeo-lab internal users:
Please see the YEOLAB_INTERNAL.md file for specific instructions on running skipper on the tscc cluster. 

# Set up
## Installation

1. **Clone the repository**  
   ```bash
   git clone https://github.com/YeoLab/skipper.git
   cd skipper
   ```

2. **Install Conda (if not already installed)**  
   The example below shows how to install Miniconda on Linux (64-bit). For detailed instructions on other systems, see the [official installation guide](https://www.anaconda.com/docs/getting-started/miniconda/install).  

   ```bash
   curl -L -O "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

3. **Create a Snakemake environment**  
   Once Conda is installed, create an environment with Snakemake version `9.12.0` to run Skipper:  

   ```bash
   conda create -n snakemake9 snakemake=9.12.0
   ```

   All other required packages and environments will be installed automatically by Snakemake and Conda the first time you run Skipper.

## Configuring Your Snakemake Profile

Snakemake profiles allow you to supply additional arguments without cluttering the command line.  
An example profile is provided at:`profiles/example_basic/config.yaml`

This profile is configured for running Skipper on a single-node machine (not recommended for most use cases; see [Running Skipper on HPCs](#Running-Skipper-on-HPCs)). The only required change is to specify a path for saving Conda environments (choose any location on your machine with sufficient storage space).

## Running Skipper on HPCs

Skipper is an end-to-end pipeline for eCLIP analysis, including:

- Preprocessing and trimming
- Alignment
- GFF partitioning
- Enriched window identification
- Fine-mapping
- Motif analysis

While Skipper can be run on powerful personal machines, it is primarily designed for **high-performance computing clusters (HPCs)**, where significant speedups are achieved by parallelizing jobs across compute nodes.

### Cluster Executor Setup

To run Skipper on HPCs, you must install a cluster executor plugin.  
The example below demonstrates installation of the **SLURM** executor (a widely used workload manager).  
Other executor options are listed in the [Snakemake plugin catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html).

```bash
conda activate snakemake9
conda install snakemake-executor-plugin-slurm=1.4.0
```

### Adjusting Your Profile

After installing the executor plugin, you must adjust your Snakemake profile.  
An example profile is provided in:

```
profiles/example_slurm/config.yaml
```

- **CONDA prefix** Do not forget to change this to a path on your cluster. 
- **Slurm Account and partition:** You must enter your own account and partition information.
- **Cluster-specific options:** Some systems require additional details. For example:  

  ```yaml
  slurm_extra: "--qos=YOUR_QOS"
  ```

  In general, if something is required in your `srun` or `sbatch` commands, it may also need to be added to `slurm_extra`.


## Minimal Example. 

This section details a small example run of skipper on a subsampled dataset. This example assumes that you are working on a linux based system with slurm set up and have already gone through all installation steps above (including adjusting the example profile). 

1. **change into your skipper directory and download the human genome from gencode**  
   ```bash
   cd annotations
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
   cd ..
   ```
2. **Edit config file**
   Open the config file in `example/Example_config.yaml` using any text editor and change every instance of /path/to/your/skipper to the **absolute** path to the skipper directory you just cloned. Also, change every instance of `/path/to/save/output` with the **absolute** path to whatever location you want to save your skipper outputs too (should be a location with lots of space, such as a scratch directory).

3. **Edit manifest file**
   Open the config file in `example/Example_config.yaml` using any text editor and change every instance of /path/to/your/skipper with the **absolute** path to the skipper directory you just cloned. 

4. **Run skipper**  
   ```bash
   unset SLURM_JOB_ID # required if running on an interactive node, which is reccomended
   
   snakemake -s Skipper.py --configfile example/Skipper_config.yaml --profile profiles/example_slurm
   ```

NOTE: The first run of skipper needs to set up all of the necessary conda environments via snakeconda and has to complete several costly steps that only need to be run for the first skipper run (e.g. parsing the GFF and generating the STAR genome index). As such, this initial skipper run will be quite slow, but subsequent runs will be much faster.

NOTE: If difficulties arrise while running this example (or any run of skipper) please see the [Troubleshooting](#Troubleshooting)) section and/or open an issue. 

# Preparing A New Skipper Run (The Config File)

Numerous resources must be entered in the `Skipper_config.yaml` file before the start of any run. These resources are split up into several different categories for ease of use:

## BASIC INPUTS
These inputs are required for all runs of skipper. 

### Required work/organizational files. 

| Resource      | Description |
| ----------- | ----------- |
| WORKDIR             | Path to save outputs to |
| TMPDIR    | Path to directory to save temporary files too (if left blank, will default to a /tmp directory inside of WORKDIR)       |
| MANIFEST            | Path to a manifest file containing information on which samples to run (Please see the "making a manifest" section)                                                      |
| TOOL_DIR    | Path to the tools directory from this repository        |

### Required Annotation Files
Each of the files in this section must already exist on your machine. Instructions for where to find/download these files for your species/cell type of interest are included in the descriptions. 

| Resource      | Description |
| ----------- | ----------- |
| GFF_source          | A short string specifying if the data came from either gencode or ensembl (options: "gencode", "ensembl") |
| GFF                 | Gzipped gene annotation to partition the transcriptome and count reads (must be from [gencode](https://www.gencodegenes.org/) or [ensembl](https://useast.ensembl.org/index.html)). |
| GENOME              | Fasta for the genome of interest (also available from gencode and ensembl) |
| ACCESSION_RANKINGS  | A ranking of gene and transcript types present in the GFF to facilitate the transcriptome partitioning  |
| BLACKLIST           | Removes windows from reproducible enriched window files. Start and end coordinates must match tiled windows exactly.  Set to None for no blacklisting    |

### Auto Generated Annotation Files.
These files can be automatically generated by Skipper. HOWEVER, it is still necessary to specify paths to these files even if they do not yet exist so that skipper knows where to save them. If these files have already been generated from other skipper runs, then specifying pre-made files will lead to significant speed-ups. 

| Input      | Description |
| ----------- | ----------- |
| PARTITION           | Gzipped BED file of windows to test (generated from GFF file) |
| FEATURE_ANNOTATIONS | Gzipped TSV file with the following columns: chrom,start,end,name,score,strand,feature_id,feature_bin,feature_type_top,feature_types,gene_name,gene_id, transcript_ids,gene_type_top,transcript_type_top,gene_types,transcript_types (generated from GFF file) |
|STAR_DIR | Directory created by STAR genomeGenerate (generated from the GFF file) |


### eCLIP Parameters. 
Each of these parameters should be edited to reflect the parameters of the specific eCLIP protocol used. 

| Setting      | Description |
| ----------- | ----------- |
| PROTOCOL            | ENCODE4 for single end, ENCODE3 for paired end |
| UMI_SIZE            | Bases to trim for deduplication (10 for current eCLIP) |
| INFORMATIVE_READ    | Which read (1 or 2) reflects the crosslink site (for Paired End runs) |
| OVERDISPERSION_MODE | Overdispersion can be estimated from multiple input replicates ("input") or multiple CLIP replicates ("clip"): "input" is recommended |
| GINI_CUTOFF: 0.9    | A filter used to remove windows with incredibly narrow peaks (skyscrapers). These skyscrapers are usually the result of PCR errors or other sequencing artifacts, and should thus be filtered out of the final result. To use a more strict filter, decrease this cutoff (e.g. 0.8). To use a more lenient filter, increase this cutoff (e.g. 0.95). To use no filter at all, just set this cutoff to any value above 1.|

## ADVANCED INPUTS
These inputs are required only if you wish to perform some of the many additional analyses available through Skipper. Removing these commands from the config file will cause Skipper to skip over these analyses.  

### Motif analysis
| Input      | Description |
| ----------- | ----------- |
| HOMER           | A boolean (True or False) specifying if you would like to run motif analysis with [Homer](http://homer.ucsd.edu/homer/motif/)|

### Meta analysis (work in progress)
| Input      | Description |
| ----------- | ----------- |
| META_ANALYSIS           | A boolean (True or False) specifying if you would like to run a meta analsis|

### Repeat analysis
| Input      | Description |
| ----------- | ----------- |
| REPEAT_TABLE | Coordinates of repetitive elements, available from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) |
| REPEAT_BED | Gzipped sorted, nonoverlapping, tab-delimited annotations of repetitive elements: chr,start,end,label,score,strand,name,class,family,proportion_gc (auto generated from repeat table) |

### Gene set enrichment analysis
| Input      | Description |
| ----------- | ----------- |
| GENE_SETS           | GMT files of gene sets for gene set enrichment calculation |
| GENE_SET_REFERENCE  | TSV of gene set name, number of windows belonging to term, and fraction of windows that lie in gene set genes |
| GENE_SET_DISTANCE   | RDS of a matrix containing jaccard index scores for all pairs of gene sets in GMT file |

# Making a manifest
The manifest file is csv a file used to direct skipper towards the raw input files for analysis. Please see the table below for a description of all required and optional columns for the table. See the example/Example_manifest.csv file for exact formatting.

NOTE: Skipper requires at least 2 replicates per sample to identify reproducible windows. 

| Column      | Description |
| ----------- | ----------- |
| Experiment       | CLIP samples will be compared against Input samples within an experiment. The same sample can be used in multiple experiments |
| Sample           | Each CLIP and Input sample will be processed separately until testing for differential binding   |
| Cells            | A place to record information on the cell sample used: this is not currently used in analysis  |
| Input_replicate  | Replicate # for the same Sample. The same Input replicate (fastq and number) can be used for multiple CLIP replicates |
| Input_adapter    | Fasta of adapter sequences for Input replicate                                                     |
| Input_fastq      | Path to Input replicate fastq (multiple files can be entered per cell to be concatenated, seperated by a space)            |
| Input_bam       | (Optional) Enter path to Input BAM file            |
| CLIP_replicate   | Replicate # for the same Sample. Distinct CLIP replicates are required |
| CLIP_adapter     | Fasta of adapter sequences for CLIP replicate                                                     |
| CLIP_fastq       | Path to CLIP replicate fastq (multiple files can be entered per cell to be concatenated, seperated by a space)            |
| CLIP_bam       | (Optional) Enter path to CLIP BAM file            |

# Skipper output 

Skipper's main outputs can be found in the WORKDIR/output folder generated by Skipper. 

Skipper produces many different outputs. The table below details some of the most important outputs that skipper can generate when all advanced inputs are specified. 

| Output      | Description |
| ----------- | ----------- |
| reproducible_enriched_windows | Table containing information/statistics of all significantly enriched windows found in all replicates |
| reproducible_enriched_re | The same as above but for repetitive regions |
| QC/multiqc/*/multiqc_report.html | Summary files that describe and visualize several quality control metrics from the initial STAR alginment of your fastq files. Very useful for confirming the integrity of your data and for observing total library size |
| homer/finemapped_results/YOUR_SAMPLE/homerResults.html |  A file that provides statistics and visuals for the top enriched binding motifs  |
| figures/reproducible_enriched_windows/*.linear.pdf | Visualization of RNA region preferences for windows called by at least two replicates   |
| figures/gene_sets        | Visualization of top enriched GO terms relative to ENCODE reproducible enriched windows   |
| figures/tsne/skipper.tsne_query.pdf       | t-SNE visualization of binding preferences releative to ENCODE RBPs   |

Skipper also creates many additional outputs and intermediate files that, while not important for basic use cases, may be useful for users
| Output      | Description |
| ----------- | ----------- |
| finemapping |
| bams |
| secondary_results/bigwigs | Bigwig files showing coverage of IP (signal) and IN (background) reads. These files can be used to observe some of the reproducible enriched windows found by skipper with a genome browser tool such as [IGV](https://igv.org/)|

WORK IN PROGRESS: Design a wiki describing all outputs. 

# Troubleshooting

1. Log files.
Skipper generates 3 types of log files. The first 2 types can be found within the `stderr/` and `stdout/` folders of your WORKDIR. These log files generally should contain all the information needed for debugging. 

However, in some cases additional information from snakemake may be necessary, in which cases users are encouraged to investigate the log files in `WORKDIR/.snakemake/slurm_logs`. These log files are organized by rules and contain additional information on the snakemake run.

2. Jobs dying with no explanation.
If you observe that many of your jobs are dying without any explanation (e.g. mostly blank files in WORKDIR/stderr, unhelpful error messages in WORKDIR/.snakemake/slurm_logs such as "Killed"), and these jobs are occuring on the same node according to WORKDIR/stdout, then it is likely that this is the result of problematic nodes on your cluster. I would reccomend taking whichever nodes were used for the failed jobs and excluding them from the analysis by adding the following lines to the slurm extra command within your profile like so:

  ```yaml
  slurm_extra: "--exclude=YOUR_NODE"
  ```

This is also the common cause of many timeout errors, as Skipper generally provides significantly more than enough time for all rules. 

3. Monitoring pipeline. 
`squeue -u $USER -o "%.18i %.10P %.20j %.10u %.2t %.10M %.6D %.20R %.80k"` will show currently active jobs (helpful to see if certain rules are getting stuck.)
