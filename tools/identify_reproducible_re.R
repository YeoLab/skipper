library(tidyverse)

# Create output directory for reproducible enriched repeat element (RE) results
dir.create("output/reproducible_enriched_re/", showWarnings = FALSE, recursive = TRUE)

# Command-line arguments:
# 1. data_directory: folder containing enriched RE TSV files
# 2. prefix: experiment prefix (used to match files)
args = commandArgs(trailingOnly=TRUE)
data_directory = args[1]
prefix = args[2]

# ---------------------------------------------
# Collect enriched RE files for this experiment
# ---------------------------------------------
# Pattern: "<prefix>.*enriched_re.tsv.gz"
# Read all non-empty tables and merge them into one tibble,
# adding "clip_replicate_label" as identifier for each replicate
enriched_re_files = list.files(path = data_directory, pattern = paste0("^", prefix, "\\..*enriched_re.tsv.gz"), full.names = TRUE)
enriched_re_data = enriched_re_files %>%
    setNames(sub("\\.enriched_windows\\.tsv.gz", "", basename(.))) %>%   # assign replicate labels from file names
    map(function(x) read_tsv(x)) %>%                                     # read each TSV
    Filter(function(x) nrow(x) > 0, .) %>%                               # drop empty files
    bind_rows(.id = "clip_replicate_label")                              # combine all replicates

# ---------------------------------------------
# Identify reproducible enriched REs
# ---------------------------------------------
# Group by all columns except replicate-specific measures
# (i.e., collapse across replicates).
# Summarize aggregate stats for each RE feature across replicates.
# Keep REs enriched (q < 0.05) in more than 1 replicate.
reproducible_enriched_re_data = enriched_re_data %>%
    group_by(across(-c(clip_replicate_label,gc_bin,baseline_l2or,input,clip,enrichment_l2or,pvalue,qvalue))) %>%
    summarize(
        input_sum = sum(input),                              # total input reads
        clip_sum = sum(clip),                                # total CLIP reads
        enrichment_n = sum(qvalue < 0.05),                   # number of replicates enriched
        enrichment_l2or_min = min(enrichment_l2or),          # min log2 enrichment
        enrichment_l2or_mean = mean(enrichment_l2or),        # mean log2 enrichment
        enrichment_l2or_max = max(enrichment_l2or),          # max log2 enrichment
        p_max = max(pvalue),                                 # worst p-value
        p_min = min(pvalue),                                 # best p-value
        q_max = max(qvalue),                                 # worst q-value
        q_min = min(qvalue)                                  # best q-value
    ) %>%
    filter(enrichment_n > 1) %>%                             # require reproducibility (>1 replicate enriched)
    arrange(p_min, desc(enrichment_l2or_mean))               # rank by strongest evidence

# ---------------------------------------------
# Save reproducible enriched REs
# ---------------------------------------------
write_tsv(reproducible_enriched_re_data, paste0("output/reproducible_enriched_re/", prefix, ".reproducible_enriched_re.tsv.gz"))
