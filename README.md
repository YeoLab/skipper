# Skipper
![Skipper cartoon](documents/logo.png)

Skip the peaks and expose RNA-binding in CLIP data

See published article in Cell Genomics: https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00085-X

## Yeo-lab internal users:
Please see the YEOLAB_INTERNAL.md file for specific instructions on running Skipper on the TSCC cluster. 

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

   Conda is used **only** to install Snakemake itself. Skipper's own dependencies do not come from Conda

4. **Make Singularity/Apptainer available**  
   Every Skipper rule runs inside one of two prebuilt container images, so you need a container runtime. On an HPC this is usually a module:

   ```bash
   module load singularitypro/3.11   # TSCC; your cluster's module name may differ
   ```

   Check it worked with `singularity --version`. Either `singularity` or `apptainer` on your `PATH` will do.

## Software environments

Every rule runs inside one of two prebuilt images published on Docker Hub:

| Image | Contents | Rules |
|---|---|---|
| `howardxu520/skipper:R_v1.0` | R 4.4.3 / Bioconductor 3.20 — tidyverse, GenomicRanges, rtracklayer, VGAM, DescTools, fgsea | The 16 `Rscript` rules |
| `howardxu520/skipper:python_v1.0` | Python 3.12 (pandas, pybedtools) plus STAR, bedtools, samtools, `bedGraphToBigWig`, fastqc, fastp, skewer, umicollapse, HOMER and MultiQC | everything from trimming and alignment through counting, coverage, bigwigs, motif calling and QC |

Snakemake pulls each image once, converts it to a `.sif`, and caches it under `apptainer-prefix`. If your compute nodes have no outbound network, pull them on a login node instead and point Skipper at the local files from your config:

```yaml
R_CONTAINER: "/abs/path/to/skipper-r.sif"
PYTHON_CONTAINER: "/abs/path/to/skipper-py.sif"
```

## Configuring Your Snakemake Profile

Snakemake profiles allow you to supply additional arguments without cluttering the command line.  
An example profile is provided at:`profiles/example_basic/config.yaml`

This profile is configured for running Skipper on a single-node machine (not recommended for most use cases; see [Running Skipper on HPCs](#Running-Skipper-on-HPCs)). Two settings need your attention:

- **`apptainer-prefix`** — where the `.sif` images are cached. Pick a location with 30 GB free that is readable from every node that will run jobs.
- **`apptainer-args`** — add a `--bind` for every directory your config points at that lives outside `WORKDIR`: `TOOL_DIR`, `GENOME`, `GFF`, `STAR_DIR`, `PARTITION`, `FEATURE_ANNOTATIONS`, `REPEAT_TABLE`, `MANIFEST` and the fastq/bam files it names. Anything not bound is simply invisible inside the container. Keep `--cleanenv`; it stops a stray `R_LIBS`, `PYTHONPATH` or `CONDA_PREFIX` in your shell from leaking in and shadowing the image's own.

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

- **`apptainer-prefix`** Do not forget to change this to a path on your cluster. It must be on shared storage that every compute node can read.
- **`apptainer-args`** Add a `--bind` for every input directory outside `WORKDIR` (see [Configuring Your Snakemake Profile](#configuring-your-snakemake-profile)).
- **Slurm account, partition:** You must enter your own account and partition information.
- **Cluster specific options:** Some systems require additional details. For example:  

  ```yaml
  slurm_extra: "--qos=YOUR_QOS"
  ```

  In general, if something is required in your `srun` or `sbatch` commands, it may also need to be added to `slurm_extra`.
- **Slurm log directory** Specifies where to save the logs from you slurm runs. Skipper generates its own log-files, but slurm logs provide a greater degree of detail. Please change this to some sensible directory on your machine. 


## Minimal Example. 

This section details a small example run of Skipper on a subsampled dataset. This example assumes that you are working on a linux based system with Slurm set up and have already gone through all installation steps above (including adjusting the example profile). 
1. **Setup an interactive node**
    While this step is technically optional, it is highly recommended to run Skipper on interactive nodes. This is especially important for your first Skipper run, when Snakemake downloads and converts the two container images. Thus, we recommend filling in the command below with your partition (-p), QOS (-q) and account (-A) information and setting up an interactive node for use with this example. Remember to `module load singularitypro/3.11` (or your cluster's equivalent) on the interactive node as well.

    ```bash
    srun -N 1 -c 1 -t 8:00:00 -p -q -A --mem 16G --pty /bin/bash
    ```

3. **Download the human genome from GENCODE**  
   ```bash
   cd /path/to/your/skipper/annotations
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
   gunzip GRCh38.primary_assembly.genome.fa.gz
   ```
4. **Edit config file**
   Open the config file in `example/Example_config.yaml` using any text editor and change every instance of `/path/to/your/skipper` to the **absolute** path to the Skipper directory you just cloned. Also, change every instance of `/path/to/save/output` with the **absolute** path to whatever location you want to save your Skipper outputs too (should be a location with lots of space, such as a scratch directory).

5. **Edit manifest file**
   Open the manifest file in `example/Example_manifest.csv` using any text editor and change every instance of /path/to/your/skipper with the **absolute** path to the Skipper directory you just cloned. 

6. **Run Skipper**  
   ```bash
   cd /path/to/your/skipper
   module load singularitypro/3.11 # or your cluster's container runtime module
   unset SLURM_JOB_ID # required if running on an interactive node. 
   snakemake -s Skipper.py --configfile example/Example_config.yaml --profile profiles/example_slurm
   ```

NOTE: On the first run, Snakemake downloads the two container images from Docker Hub and converts them to `.sif` files under your `apptainer-prefix`. That run also has to complete several costly steps that only need to be run once (e.g. parsing the GFF and generating the STAR genome index). As such, this initial Skipper run will be quite slow, but subsequent runs reuse the cached images and start immediately.

NOTE: If difficulties arrise while running this example (or any run of Skipper) please see the [Troubleshooting](#Troubleshooting) section and/or open an issue. 

# Preparing A New Skipper Run 

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

### Running the GFF filtration

Skipper works with GFF files from either GENCODE or Ensembl. This section provides a brief example of running the filtration pipeline on the `gencode.v49.basic.annotation.gff3.gz` file available from the [GENCODE website](https://www.gencodegenes.org/human/).

The GFF utilities are not Snakemake rules — you run them by hand, once, before
your first Skipper run — but they need no Conda environment either. Both
Skipper images carry their dependencies, so run them with `apptainer exec`.

1. **Point at the images**
    ```bash
    module load singularitypro/3.11   # or your cluster's container runtime module
    cd path/to/your/skipper/gff_utils

    # Either the .sif files you built, or ones pulled from Docker Hub:
    #   apptainer pull skipper-r.sif  docker://howardxu520/skipper:R_v1.0
    #   apptainer pull skipper-py.sif docker://howardxu520/skipper:python_v1.0
    SKIPPER_R=/abs/path/to/skipper-r.sif
    SKIPPER_PY=/abs/path/to/skipper-py.sif
    ```

2. **Run the filtration pipeline**
    ```bash
    apptainer exec --cleanenv --bind "$PWD" "$SKIPPER_R" \
      Rscript "path/to/your/gff_transcript_quality_filter.R" \
        "gencode" \
        "path/to/your/gencode.v49.basic.annotation.gff3.gz" \
        "./gencode.v49.filtered.annotation.gff3.gz"
    ```

    Add a `--bind` for the directory holding your GFF if it lives outside the
    current one; only paths you bind are visible inside the container.

And that's it. To use an Ensembl GFF instead of a GENCODE GFF, simply replace the `"gencode"` argument supplied to the R script with `"ensembl"`.

#### Optional: Filtering by expression levels

Some users may wish to further filter their GFF files using known expression levels for their cell type. This can increase the power of downstream analyses, as Skipper will not evaluate genes or transcripts that are known to be unexpressed in the cell type of interest.

Skipper accepts three sources of expression data:

- [GTEx](https://gtexportal.org/home/) (human tissue types)
- [CCLE](https://sites.broadinstitute.org/ccle) (common cell lines)
- [Salmon](https://combine-lab.github.io/salmon/) outputs (custom/user-generated)

The example below shows how to filter for genes with a TPM value greater than 1 using GTEx data from Breast Mammary Tissue samples, available from the [GTEx website](https://www.gtexportal.org/home/downloads/adult-gtex%23qtl). Note that this filtration is performed on a GFF that has already gone through the quality-filtering step described above.

3. **Filter a GFF by expression level**

    This one runs in the Python image rather than the R image, since it is a
    `pyranges` script:

    ```bash
    apptainer exec --cleanenv --bind "$PWD" --bind "${Skipper_dir}" "$SKIPPER_PY" \
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
|NORMALIZATION_MODE | Use Skipper's new pseudocount based normalization mode ("new") or its original expected ratio adjustment normalization ("classic"). "new" is reccomended. Classic is maintained only to ensure reproducibility for older analyses. |
|THRESHOLD_MIN | The minimum threshold strength used by Flipper. Setting this to N will prevent Skipper from identifying windows with less than N total reads across the input and IP fractions as significantly enriched. |
|DEFINITION_OF_REPRODUCIBILITY | Minimum number of replicates with significant enrichment to be considered "reproducible"*. Use "ALL" to automatically select the total number of replicates for the sample as the cutoff (the default behavior), otherwise enter an integer.|

*Skipper sticks to a traditional definition of "reproducible" where any window found as significant in all replicates is considered reproducible. However, some users may want to introduce a less strict cutoff when working with higher replicates, such as any window found in 2/3 replicates.

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

2. Container problems.

- `The apptainer or singularity command has to be available...`, you need to make sure the singularity is setup on your system.
- `No such file or directory` for an input that plainly exists, the path is almost certainly not bound into the container. Add a `--bind` for it to `apptainer-args` in your profile. This is the single most common failure when migrating a working config: only `WORKDIR` is bound automatically.
- The image download fails on compute nodes, they likely have no outbound network. Pull the images on a login node and set `R_CONTAINER` / `PYTHON_CONTAINER` in your Skipper config to the resulting `.sif` paths.
- An R or Python package appears to be the wrong version, a host `R_LIBS`, `R_LIBS_USER` or `PYTHONPATH` is leaking in. Make sure `--cleanenv` is present in `apptainer-args`.

3. Jobs dying with no explanation.
If you observe that many of your jobs are dying without any explanation (e.g. mostly blank files in WORKDIR/stderr, unhelpful error messages in WORKDIR/.snakemake/slurm_logs such as "Killed"), and these jobs are occuring on the same node according to WORKDIR/stdout, then it is likely that this is the result of problematic nodes on your cluster. I would reccomend taking whichever nodes were used for the failed jobs and excluding them from the analysis by adding the following lines to the slurm extra command within your profile like so:

  ```yaml
  slurm_extra: "--exclude=YOUR_NODE"
  ```

This is also the common cause of many timeout errors, as Skipper generally provides significantly more than enough time for all rules. 

4. Monitoring pipeline. 
`squeue -u $USER -o "%.18i %.10P %.20j %.10u %.2t %.10M %.6D %.20R %.80k"` will show currently active jobs (helpful to see if certain rules are getting stuck.)
