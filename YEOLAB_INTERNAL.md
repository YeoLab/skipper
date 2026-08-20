# YEO-LAB internal example 
Hello. This is a short example for running Skipper as a member of the Yeo-lab partition on TSCC.

## Load up an interactive node.

```
srun -N 1 -c 1 -t 8:00:00 -p gold -q hcg-csd792 -A csd792 --mem 16G --pty /bin/bash
```

## Create a folder in scratch to save the output.
Run the command below (Replacing YOUR_USERNAME with your TSCC username) to create a folder to save the output of the run. 

```
mkdir /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper2_test
cd /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper2_test
```

## Load up the skipper module.

```
module load skipper
```

**IMPORTANT**: The only modules that should be loaded are the default modules (the ones that automatically load whenever you go on TSCC) and Skipper (Loading Skipper also automatically loads the singularity module). No other modules should be loaded, as this may confuse snakemake as to which environment to run the code in, causing an error. You can check what modules you have loaded using ```module list ```

## Copy the example config file from the Skipper module. 

Copy the config file from the skipper module folder to the scratch directory you just made. 

```
cp $SKIPPER_HOME/bin/skipper/example/yeo_lab_internal_example_config.yaml ./yeo_lab_internal_example_config.yaml
```

##  Adjust the yeo loab internal example config file. 
Adjust the WOKDIR input in your copy of `yeo_lab_internal_example_config.yaml` so that it points to the scratch directory you made in the first step. The config file can be adjsuted using any text editor you are comfortable with (vim, nano, jupyternotebooks, etc). 

No other changes to the config are necessary. 

## Run skipper. 

Now, simply replace YOUR_USERNAME in the command below again and run the following commands. 

```
unset SLURM_JOB_ID
   
snakemake -s $SKIPPER_HOME/bin/skipper/Skipper.py --configfile /tscc/lustre/ddn/scratch/YOUR_USERNAME/skipper2_test/yeo_lab_internal_example_config.yaml --profile $SKIPPER_HOME/bin/skipper/profiles/tscc2_snakemake9
```

**NOTE:** If problems occur in this initial run (or any run), please check the troubleshooting section below (especially point 2). If problems persist, please open up a github issue at [https://github.com/YeoLab/skipper](https://github.com/YeoLab/skipper)

# YEO-LAB internal notes on preparing a new run. 
Running skipper on TSCC with a new eCLIP dataset generally only requires 2 things, a config file and a manifest. 

Below I have copied and modified several sections from the main readme on making a configfile and manifest, along with a short section on troubleshooting. I have added Yeo-lab internal user specific notes where appropriate under sections labelled **YEOLAB INTERNAL USER NOTE**.

# Preparing A New Skipper Run (The Config File)

## Filtering/pre-processing GFF files

Skipper uses GFF files both to provide users with metadata on the types of genomic features an RBP is binding (e.g., introns, exons, UTRs, etc.) and to guide the creation of genomic windows to test (e.g., Skipper attempts to avoid creating windows that cross intron/exon boundaries).

As such, it is incredibly important to ensure that these GFF files have been properly cleaned and filtered before supplying them to Skipper. GFF files, by design, include many transcripts from alternative isoforms. While these transcripts are biologically meaningful, it is necessary to select one "best" transcript per gene to ensure that Skipper's results are interpretable and that its window selection is optimal.

The GFF filtration method provided with Skipper consists of two parts. First, transcripts are filtered using the "tag" column of the GFF file, where we attempt to select the best transcript for each gene based on transcript-quality annotations provided by GENCODE or Ensembl. For example, GENCODE transcripts are ranked according to the following hierarchy:

| Rank | Annotation |
|------|------------|
| 1 | MANE_Select |
| 2 | Ensembl_canonical |
| 3 | GENCODE_Primary |
| 4 | APPRIS principal |
| 5 | CCDS |
| 6 | GENCODE basic |

When multiple transcripts receive the same highest-ranking annotation for a gene, additional transcript-quality annotations are used as tie-breakers. If a tie still remains after all available transcript-quality information has been considered, all tied transcripts are retained.

Second, Skipper resolves overlapping feature annotations. Even after transcript filtering, some genomic regions may still be assigned to multiple feature types (e.g., a region may be annotated as both an exon and an intron due to differences between retained transcript isoforms). In these cases, Skipper assigns a primary feature type using the accession type ranking file provided to the program.

It is important to note that although we feel this is the best GFF filtering strategy for most use cases, there are scenarios where such aggressive GFF filtering may be counterproductive. For example, a user may be interested in RBP binding to rare isoforms, or a particular treatment may have substantially altered the isoform distribution within the cells being studied. In such cases, we recommend that users generate their own custom GFF files.

**YEOLAB INTERNAL USER NOTE:**

Yeo-lab members on TSCC have access to several pre-built annotation resources, including GFF, PARTITION, FEATURE_ANNOTATION, and STAR_DIR files. These are available in `/tscc/projects/ps-yeolab4/software/skipper/2.0.0/bin/skipper/annotations`

Using these whenever possible will lead to significant speedups. **DO NOT USE ANY ANNOTATION FILES FROM OLDER SKIPPER RUNS!!!!** these files were generated before the implementation of the GFF file filtration (removes problematic transcripts) and can lead to erroneous results. 

### Running the GFF filtration

Skipper works with GFF files from either GENCODE or Ensembl. This section provides a brief example of running the filtration pipeline on the `gencode.v49.basic.annotation.gff3.gz` file available from the [GENCODE website](https://www.gencodegenes.org/human/).

1. **Create and activate a new conda environment from the provided YAML file**
    ```bash
    cd path/to/your/skipper/gff_utils
    conda env create -f gff_utils.yaml
    conda activate gff_utils_env
    ```

2. **Run the filtration pipeline**
    ```bash
    Rscript "path/to/your/gff_transcript_quality_filter.R" \
        "gencode" \
        "path/to/your/gencode.v49.basic.annotation.gff3.gz" \
        "./gencode.v49.filtered.annotation.gff3.gz"
    ```

And that's it. To use an Ensembl GFF instead of a GENCODE GFF, simply replace the `"gencode"` argument supplied to the R script with `"ensembl"`.

#### Optional: Filtering by expression levels

Some users may wish to further filter their GFF files using known expression levels for their cell type. This can increase the power of downstream analyses, as Skipper will not evaluate genes or transcripts that are known to be unexpressed in the cell type of interest.

Skipper accepts three sources of expression data:

- [GTEx](https://gtexportal.org/home/) (human tissue types)
- [CCLE](https://sites.broadinstitute.org/ccle) (common cell lines)
- [Salmon](https://combine-lab.github.io/salmon/) outputs (custom/user-generated)

The example below shows how to filter for genes with a TPM value greater than 1 using GTEx data from Breast Mammary Tissue samples, available from the [GTEx website](https://www.gtexportal.org/home/downloads/adult-gtex%23qtl). Note that this filtration is performed on a GFF that has already gone through the quality-filtering step described above.

3. **Filter a GFF by expression level**
    ```bash
    python "${Skipper_dir}/gff_utils/gff_expression_filter.py" \
        -a "path/to/your/gencode.v49.filtered.annotation.gff3.gz" \
        -t 1 \
        -s "GTEx" \
        -c "Breast_Mammary_Tissue" \
        -q "path/to/your/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct" \
        -o "./gencode.GTEX_test.gz"
    ```

## The Config File

Numerous resources must be entered in the `Skipper_config.yaml` file before the start of any run. These resources are split up into several different categories for ease of use:

### Basic inputs
These inputs are required for all runs of Skipper. 

#### Required work/organizational files

| Resource      | Description |
| ----------- | ----------- |
| WORKDIR             | Path to save outputs to |
| TMPDIR    | Path to directory to save temporary files too (if left blank, will default to a /tmp directory inside of WORKDIR)       |
| MANIFEST            | Path to a manifest file containing information on which samples to run (Please see the "making a manifest" section)                                                      |
| TOOL_DIR    | Path to the tools directory from this repository        |

#### Required Annotation Files
Each of the files in this section must already exist on your machine. Instructions for where to find/download these files for your species/cell type of interest are included in the descriptions. 

| Resource      | Description |
| ----------- | ----------- |
| GFF                 | Gzipped gene annotation to partition the transcriptome and count reads (must be from [GENCODE](https://www.gencodegenes.org/) or [Ensembl](https://useast.ensembl.org/index.html)). |
| GENOME              | FASTA for the genome of interest (also available from GENCODE and Ensembl) |
| ACCESSION_RANKINGS  | A ranking of gene and transcript types present in the GFF to facilitate the transcriptome partitioning  |
| BLACKLIST           | Removes windows from reproducible enriched window files. Start and end coordinates must match tiled windows exactly.  Leave blank for no blacklist    |

NOTE: All gene and transcript types present in the GFF file must also be present in the ACCESSION_RANKINGS file.

**YEOLAB INTERNAL USER NOTE**
The parameters GENOME, ACCESSION_RANKINGS, BLACKLIST, REPEAT_TABLE, REPEAT_BED, GENE_SETS, GENE_SET_REFERENCE, and GENE_SET_DISTANCE can typically remain unchanged from the settings in `yeo_lab_internal_example_config.yaml`.

#### Auto Generated Annotation Files
These files can be automatically generated by Skipper. HOWEVER, it is still necessary to specify paths to these files even if they do not yet exist so that Skipper knows where to save them. If these files have already been generated from other Skipper runs, then specifying pre-made files will lead to significant speed-ups. 

| Input      | Description |
| ----------- | ----------- |
| PARTITION           | Gzipped BED file of windows to test (generated from GFF file) |
| FEATURE_ANNOTATIONS | Gzipped TSV file with the following columns: chrom,start,end,name,score,strand,feature_id,feature_bin,feature_type_top,feature_types,gene_name,gene_id, transcript_ids,gene_type_top,transcript_type_top,gene_types,transcript_types (generated from GFF file) |
|STAR_DIR | Directory created by STAR genomeGenerate (generated from the GFF file) |


#### eCLIP Parameters
Each of these parameters should be edited to reflect the parameters of the specific eCLIP protocol used. 

| Setting      | Description |
| ----------- | ----------- |
| PROTOCOL            | "ENCODE4" for single end, "ENCODE3" for paired end (NOTE: Only used when running Skipper with FASTQs. will be ignored if running Skipper using BAM files)|
| UMI_SIZE            | Bases to trim for deduplication (10 for current eCLIP) (NOTE: Only used when running Skipper with FASTQs. will be ignored if running Skipper using BAM files) |
| INFORMATIVE_READ    | Which read (1 or 2) reflects the crosslink site (for Paired End runs) |
| OVERDISPERSION_MODE | Overdispersion can be estimated from multiple input replicates ("input") or multiple CLIP replicates ("clip"): "input" is recommended |
| GINI_CUTOFF   | A filter used to remove windows with incredibly narrow peaks (skyscrapers). These skyscrapers are usually the result of PCR errors or other sequencing artifacts, and should thus be filtered out of the final result. To use a more strict filter, decrease this cutoff (e.g. 0.8). To use a more lenient filter, increase this cutoff (e.g. 0.95). To use no filter at all, just set this cutoff to any value above 1.|
|NORMALIZATION_MODE | Use Skipper's new pseudocount based normalization mode ("new") or its original expected ratio adjustment normalization ("classic"). "new" is reccomended. Classic is maintained only to ensure reproducibility for older analyses. 
|THRESHOLD_MIN | The minimum threshold strength used by Flipper. Setting this to N will prevent Skipper from identifying windows with less than N total reads across the input and IP fractions as significantly enriched. 

**YEOLAB INTERNAL USER NOTE**
Most eCLIP parameters (e.g., protocol, UMI_SIZE, GINI_CUTOFF) can remain at their default values from `yeo_lab_internal_example_config.yaml` for the majority of datasets.

### Advanced inputs
These inputs are required only if you wish to perform some of the many additional analyses available through Skipper. Removing these commands from the config file will cause Skipper to skip over these analyses.  

#### Motif analysis
| Input      | Description |
| ----------- | ----------- |
| HOMER           | A boolean (True or False) specifying if you would like to run motif analysis with [Homer](http://homer.ucsd.edu/homer/motif/)|

#### Repeat analysis
| Input      | Description |
| ----------- | ----------- |
| REPEAT_TABLE | Coordinates of repetitive elements, available from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) |
| REPEAT_BED | Gzipped sorted, nonoverlapping, tab-delimited annotations of repetitive elements: chr,start,end,label,score,strand,name,class,family,proportion_gc (auto generated from repeat table) |

#### Gene set enrichment analysis
| Input      | Description |
| ----------- | ----------- |
| GENE_SETS           | GMT files of gene sets for gene set enrichment calculation |
| GENE_SET_REFERENCE  | TSV of gene set name, number of windows belonging to term, and fraction of windows that lie in gene set genes |
| GENE_SET_DISTANCE   | RDS of a matrix containing jaccard index scores for all pairs of gene sets in GMT file |

## Making a manifest
The manifest file is a csv file used to direct Skipper to the input files for analysis. The table below lists required annotation information that MUST be included in the every manifest.

NOTE: Skipper requires at least 2 replicates per sample to identify reproducible windows. 
NOTE: See the example/Example_manifest.csv file for example formatting.

| Column      | Description |
| ----------- | ----------- |
| Experiment       | CLIP samples will be compared against Input samples within an experiment. The same sample can be used in multiple experiments |
| Sample           | Each CLIP and Input sample will be processed separately until testing for differential binding   |
| Cells            | A place to record information on the cell sample used: this is not currently used in analysis  |
| Input_replicate  | Replicate number for the same Sample. The same Input replicate (FASTQ and number) can be used for multiple CLIP replicates |
| CLIP_replicate   | Replicate number for the same Sample. Distinct CLIP replicates are required |

The remaining required inputs for the manifest depend on if you are running Skipper with FASTQ files or BAM files. 

## Running Skipper with FASTQs

Add the following columns to your manifest
| Column      | Description |
| ----------- | ----------- |
| Input_adapter    | Path to FASTA file containing adapter sequences for Input replicate                                                     |
| Input_fastq      | Path to Input replicate FASTQ (multiple files can be entered per cell to be concatenated, seperated by a space)            |
| CLIP_adapter     | Path to FASTA file containing adapter sequences for CLIP replicate                                                     |
| CLIP_fastq       | Path to CLIP replicate FASTQ (multiple files can be entered per cell to be concatenated, seperated by a space)            |

## Running Skipper with BAMs

Add the following columns to your manifest
| Column      | Description |
| ----------- | ----------- |
| Input_bam       | Path to Input BAM file            |
| CLIP_bam       |  Path to CLIP BAM file            |

NOTE: When runnign Skipper with BAM files, it is important to ensure that the GFF file used when creating the BAM files is identical to the one given as input to Skipper. 
NOTE: When using BAM files, outputs related to Skipper’s pre-processing steps (e.g., quality control analyses) will not be generated.

# Skipper output 

Skipper's main outputs can be found in the WORKDIR/output folder generated by Skipper. 

Skipper produces many different outputs. The table below details some of the most important outputs that Skipper can generate when all advanced inputs are specified. 

| Output      | Description |
| ----------- | ----------- |
| reproducible_enriched_windows | Table containing information/statistics of all significantly enriched windows found in all replicates |
| reproducible_enriched_re | The same as above but for repetitive regions |
| QC/multiqc/*/multiqc_report.html | Summary files that describe and visualize several quality control metrics from the initial STAR alginment of your FASTQ files. Very useful for confirming the integrity of your data and for observing total library size |
| homer/finemapped_results/YOUR_SAMPLE/homerResults.html |  A file that provides statistics and visuals for the top enriched binding motifs  |
| figures/reproducible_enriched_windows/*.linear.pdf | Visualization of RNA region preferences for windows called by at least two replicates   |
| figures/gene_sets        | Visualization of top enriched GO terms relative to ENCODE reproducible enriched windows   |
| figures/tsne/skipper.tsne_query.pdf       | t-SNE visualization of binding preferences releative to ENCODE RBPs   |

Skipper also creates many additional outputs and intermediate files that, while not needed for most use cases, may be useful for some users. 
| Output      | Description |
| ----------- | ----------- |
| secondary_results/bams | BAM files generated using Skipper's native pre-processing pipeline. |
| secondary_results/bigwigs | Bigwig files showing coverage of IP (signal) and IN (background) reads. These files can be used to observe some of the reproducible enriched windows found by Skipper with a genome browser tool such as [IGV](https://igv.org/)|

# Troubleshooting

1. Log files.
Skipper generates 3 types of log files. The first 2 types can be found within the `stderr/` and `stdout/` folders of your WORKDIR. These log files generally should contain all the information needed for debugging. 

However, in some cases additional information from snakemake may be necessary, in which cases users are encouraged to investigate the log files in `WORKDIR/.snakemake/slurm_logs`. These log files are organized by rules and contain additional information on the snakemake run.

2. snakeconda installation.
When running Skipper for the first time, it is not uncommon to run into `CreateCondaEnvironmentException:` errors. While their are a variety of factors that can contribute to these errors (make sure your conda is updated) we have found that these errors can usually be attributed to snakeconda running out of ram when setting up some of the heavier environemnts. Switching to an interactive node with more ram usually solves the issue. 

3. Jobs dying with no explanation.
If you observe that many of your jobs are dying without any explanation (e.g. mostly blank files in WORKDIR/stderr, unhelpful error messages in WORKDIR/.snakemake/slurm_logs such as "Killed"), and these jobs are occuring on the same node according to WORKDIR/stdout, then it is likely that this is the result of problematic nodes on your cluster. I would reccomend taking whichever nodes were used for the failed jobs and excluding them from the analysis by adding the following lines to the slurm extra command within your profile like so:

  ```yaml
  slurm_extra: "--exclude=YOUR_NODE"
  ```

This is also the common cause of many timeout errors, as Skipper generally provides significantly more than enough time for all rules. 

4. Monitoring pipeline. 
`squeue -u $USER -o "%.18i %.10P %.20j %.10u %.2t %.10M %.6D %.20R %.80k"` will show currently active jobs (helpful to see if certain rules are getting stuck.)