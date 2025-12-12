# YEO-LAB internal example 
Hello. This is a short example for running skipper as a member of the Yeo-lab partition on TSCC. Before attempting this example, please log onto TSCC and change directories to your scratch director

## Create a folder in scratch to save the output.
After logging onto TSCC, simply run the command below (Replacing YOUR_USERNAME with your TSCC username) to create a folder to save the output of the 

```
mkdir /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper100_test
cd /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper100_test
```

## Load up the skipper module.

```
module load skipper
```

## Download the example files from github. 

Copy the config file from the skipper module folder to somewhere on your TSCC (e.g. a scratch direcotry). 

```
cp $SKIPPER_HOME/bin/skipper/example/yeo_lab_internal_example_config.yaml ./yeo_lab_internal_example_config.yaml
```

##  Adjust the yeo loab internal example config file. 
Adjust the WOKDIR input in your copy of `yeo_lab_internal_example_config.yaml` so that it points to the scratch directory you made in the first step. The config file can be adjsuted using any text editor you are comfortable with (vim, nano, jupyternotebooks, etc). 

No other changes to the config are necessary. 

## Run skipper. 

Now, simply replace YOUR_USERNAME in the command below again and run the following commands. 

```
srun -N 1 -c 1 -t 4:00:00 -p gold -q hcg-csd792 -A csd792 --mem 4G --pty /bin/bash

unset SLURM_JOB_ID
   
snakemake -s $SKIPPER_HOME/bin/skipper/Skipper.py --configfile /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper100_test/yeo_lab_internal_example_config.yaml --profile $SKIPPER_HOME/bin/skipper/profiles/tscc2_snakemake9
```

**NOTE:** If problems occur in this initial run (or any run), please check the troubleshooting section below (especially point 2). If problems persist, please open up a github issue at [https://github.com/YeoLab/skipper](https://github.com/YeoLab/skipper)

# YEO-LAB internal notes on preparing a new run. 
Running skipper on TSCC with a new eCLIP dataset generally only requires 2 things, a config file and a manifest. 

Below I have copied and modified several sections from the main readme on making a configfile and manifest, along with a short section on troubleshooting. I have added Yeo-lab internal user specific notes where appropriate under sections labelled **YEOLAB INTERNAL USER NOTE**.

# Preparing A New Skipper Run (The Config File)

Numerous resources must be entered in the `Skipper_config.yaml` file before the start of any run. These resources are split up into several different categories for ease of use:

**YEOLAB INTERNAL USER NOTE**
1. Yeo-lab members on TSCC have access to several pre-built annotation resources, including GFF, PARTITION, FEATURE_ANNOTATION, and STAR_DIR files. These are available in `/tscc/projects/ps-yeolab4/software/skipper/1.100.0/bin/skipper/annotations`

   
    Using these whenever possible will lead to significant speedups. **DO NOT USE ANY ANNOTATION FILES FROM OLDER SKIPPER RUNS!!!!** these files were generated before the implementation of the GFF file filtration (removes problematic transcripts) and can lead to erroneous results. 

3. The parameters GENOME, ACCESSION_RANKINGS, BLACKLIST, REPEAT_TABLE, REPEAT_BED, GENE_SETS, GENE_SET_REFERENCE, and GENE_SET_DISTANCE can typically remain unchanged from the settings in `yeo_lab_internal_example_config.yaml`.
4. Most eCLIP parameters (e.g., protocol, UMI_SIZE, GINI_CUTOFF) can remain at their default values from `yeo_lab_internal_example_config.yaml` for the majority of datasets.

## BASIC INPUTS
These inputs are required for all runs of skipper. 

### Required work/organizational files. 

| Resource      | Description |
| ----------- | ----------- |
| WORKDIR   | Path to save outputs to |
| TMPDIR    | Path to directory to save temporary files too (if left blank, will default to a /tmp directory inside of WORKDIR) |
| MANIFEST  | Path to a manifest file containing information on which samples to run (Please see the "making a manifest" section) |
| TOOL_DIR  | Path to the tools directory from this repository |

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
|STAR_DIR             | Directory created by STAR genomeGenerate (generated from the GFF file) |


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
| HOMER | A boolean (True or False) specifying if you would like to run motif analysis with [Homer](http://homer.ucsd.edu/homer/motif/)|

### Meta analysis (work in progress)
| Input      | Description |
| ----------- | ----------- |
| META_ANALYSIS | A boolean (True or False) specifying if you would like to run a meta analsis|

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
The manifest file is csv a file used to direct skipper towards the raw input files for analysis. Please see the table below for a description of all required and optional columns for the table. See the `$SKIPPER_HOME/bin/skipper/example/yeo_lab_internal_example_manifest.csv` file for an example of the exact formatting (this is the same manifest used in our example).

NOTE: Skipper requires at least 2 replicates per sample to identify reproducible windows. 

| Column      | Description |
| ----------- | ----------- |
| Experiment       | CLIP samples will be compared against Input samples within an experiment. The same sample can be used in multiple experiments |
| Sample           | Each CLIP and Input sample will be processed separately until testing for differential binding   |
| Cells            | A place to record information on the cell sample used: this is not currently used in analysis  |
| Input_replicate  | Replicate # for the same Sample. The same Input replicate (fastq and number) can be used for multiple CLIP replicates |
| Input_adapter    | Fasta of adapter sequences for Input replicate |
| Input_fastq      | Path to Input replicate fastq (multiple files can be entered per cell to be concatenated, seperated by a space) |
| Input_bam        | (Optional) Enter path to Input BAM file |
| CLIP_replicate   | Replicate # for the same Sample. Distinct CLIP replicates are required |
| CLIP_adapter     | Fasta of adapter sequences for CLIP replicate |
| CLIP_fastq       | Path to CLIP replicate fastq (multiple files can be entered per cell to be concatenated, seperated by a space) |
| CLIP_bam         | (Optional) Enter path to CLIP BAM file |

# Troubleshooting

1. Log files.
Skipper generates 3 types of log files. The first 2 types can be found within the `stderr/` and `stdout/` folders of your WORKDIR. These log files generally should contain all the information needed for debugging. That being said, in some cases additional information from snakemake may be necessary, in which case users are encouraged to investigate the log files in `WORKDIR/.snakemake/slurm_logs`. These log files are organized by rules and contain additional information on the snakemake run.

2. Jobs dying with no explanation.
If you observe that many of your jobs are dying without any explanation (e.g. mostly blank files in WORKDIR/stderr, unhelpful error messages in WORKDIR/.snakemake/slurm_logs such as "Killed"), and these jobs are occuring on the same node according to WORKDIR/stdout, then it is likely that this is the result of problematic nodes on your cluster, and I would reccomend the following steps:

    - copy the profile snakemake uses to communicate with slurm to somewhere on TSCC:
      ```
      cp $SKIPPER_HOME/bin/skipper/profiles/tscc2_snakemake9/config.yaml ANYWHERE/ON/TSCC/profiles/tscc2_snakemake9/config.yaml
      ```
    - take whichever nodes were used for the failed jobs and excluding them from the analysis by adding them to the exclude section of your copied profile like so:
      ```
      slurm_extra: "--qos=hcg-csd792 --exclude=tscc-1-18, YOUR_PROBLEM_NODE"
      ```
        
    - Replace the profile in the run with your adjusted profile:

        ```
      snakemake -s $SKIPPER_HOME/bin/skipper/Skipper.py --configfile /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper100_test/yeo_lab_internal_example_config.yaml --profile ANYWHERE/ON/TSCC/profiles/tscc2_snakemake9
        ```
    
This is also the common cause of many timeout errors, as Skipper generally provides significantly more than enough time for all rules. 

3. Monitoring pipeline. 
`squeue -u $USER -o "%.18i %.10P %.20j %.10u %.2t %.10M %.6D %.20R %.80k"` will show currently active jobs (helpful to see if certain rules are getting stuck.)
