library(tidyverse)

# Create output directory for reproducible enriched repeat element (RE) results.
dir.create("output/reproducible_enriched_re/", showWarnings = FALSE, recursive = TRUE)

# Also create a directory for sentinel plots for reproducible enriched REs.
dir.create("output/figures/secondary_figures/reproducible_enriched_re/",
           showWarnings = FALSE, recursive = TRUE)

# Command-line arguments:
# 1. data_directory: folder containing enriched RE TSV files
# 2. prefix: experiment prefix (used to match files)
args = commandArgs(trailingOnly=TRUE)
data_directory = args[1]
prefix = args[2]

# ---------------------------------------------
# Collect enriched RE files for this experiment.
# ---------------------------------------------
enriched_re_files = list.files(
    path   = data_directory,
    pattern = paste0("^", prefix, "\\..*enriched_re.tsv.gz"),
    full.names = TRUE
)

# Function to detect a dummy enriched_re file.
is_dummy_enriched_re <- function(df) {
    if (nrow(df) == 0) return(TRUE)

    key_numeric_cols = intersect(
        c("input", "clip", "enrichment_l2or", "pvalue", "qvalue"),
        colnames(df)
    )

    # If all key numeric cols are NA for *all* rows → dummy.
    all(sapply(df[key_numeric_cols], function(x) all(is.na(x))))
}

# Load tables, remove dummy ones.
enriched_re_list =
    enriched_re_files %>%
    setNames(sub("\\.enriched_windows\\.tsv.gz", "", basename(.))) %>%
    map(read_tsv) %>%
    discard(is_dummy_enriched_re)

# If all files were dummy → result is empty.
if (length(enriched_re_list) == 0) {
    enriched_re_data = tibble()
} else {
    enriched_re_data = bind_rows(enriched_re_list, .id = "clip_replicate_label")
}

# ---------------------------------------------
# Identify reproducible enriched REs.
# ---------------------------------------------

if (nrow(enriched_re_data) > 0) {
    reproducible_enriched_re_data = enriched_re_data %>%
        group_by(across(-c(clip_replicate_label,
                           gc_bin,
                           baseline_l2or,
                           input,
                           clip,
                           enrichment_l2or,
                           pvalue,
                           qvalue))) %>%
        summarize(
            input_sum = sum(input, na.rm = TRUE),              # Total input reads.
            clip_sum = sum(clip, na.rm = TRUE),                # Total CLIP reads.
            enrichment_n = sum(qvalue < 0.05, na.rm = TRUE),   # Number of replicates enriched.
            enrichment_l2or_min = min(enrichment_l2or, na.rm = TRUE),
            enrichment_l2or_mean = mean(enrichment_l2or, na.rm = TRUE),
            enrichment_l2or_max = max(enrichment_l2or, na.rm = TRUE),
            p_max = max(pvalue, na.rm = TRUE),
            p_min = min(pvalue, na.rm = TRUE),
            q_max = max(qvalue, na.rm = TRUE),
            q_min = min(qvalue, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        # If everything was NA, the *_min/max/mean can become Inf/-Inf; coerce to NA.
        mutate(
            across(c(enrichment_l2or_min,
                     enrichment_l2or_mean,
                     enrichment_l2or_max,
                     p_max,
                     p_min,
                     q_max,
                     q_min),
                   ~ ifelse(is.infinite(.), NA_real_, .))
        ) %>%
        filter(enrichment_n > 1) %>%                             # Require reproducibility (>1 replicate enriched).
        arrange(p_min, desc(enrichment_l2or_mean))               # Rank by strongest evidence.

} else {
    # No enriched_re_data at all; start with an empty result.
    reproducible_enriched_re_data = tibble()
}

# ---------------------------------------------
# If no reproducible REs, create a dummy table.
# ---------------------------------------------
if (nrow(reproducible_enriched_re_data) < 1) {
    # Derive grouping columns from enriched_re_data, if available.
    base_group_cols = setdiff(
        names(enriched_re_data),
        c("clip_replicate_label",
          "gc_bin",
          "baseline_l2or",
          "input",
          "clip",
          "enrichment_l2or",
          "pvalue",
          "qvalue")
    )

    if (nrow(enriched_re_data) > 0 && length(base_group_cols) > 0) {
        # Use the first row to preserve types for grouping columns.
        reproducible_enriched_re_data = enriched_re_data %>%
            slice(1) %>%
            select(all_of(base_group_cols)) %>%
            mutate(
                input_sum = NA_real_,
                clip_sum = NA_real_,
                enrichment_n = NA_integer_,
                enrichment_l2or_min = NA_real_,
                enrichment_l2or_mean = NA_real_,
                enrichment_l2or_max = NA_real_,
                p_max = NA_real_,
                p_min = NA_real_,
                q_max = NA_real_,
                q_min = NA_real_
            )
    } else {
        # Full fallback if we truly have no structure to borrow.
        reproducible_enriched_re_data = tibble(
            input_sum = NA_real_,
            clip_sum = NA_real_,
            enrichment_n = NA_integer_,
            enrichment_l2or_min = NA_real_,
            enrichment_l2or_mean = NA_real_,
            enrichment_l2or_max = NA_real_,
            p_max = NA_real_,
            p_min = NA_real_,
            q_max = NA_real_,
            q_min = NA_real_
        )
    }

    # Create a sentinel plot indicating that there were no reproducible enriched REs.
    sentinel_pdf = file.path(
        "output/figures/secondary_figures/reproducible_enriched_re/",
        paste0(prefix, ".reproducible_enriched_re_dummy.pdf")
    )
    pdf(sentinel_pdf, height = 1.8, width = 2.8)
    par(mar = c(0,0,0,0))      # Critical fix: remove default margins.
    plot.new()
    text(0.5, 0.5,
         labels = "No reproducible enriched repeat elements.\nDummy output generated.",
         cex = 0.7)
    dev.off()
}

# ---------------------------------------------
# Save reproducible enriched REs (real or dummy).
# ---------------------------------------------
write_tsv(
    reproducible_enriched_re_data,
    paste0("output/reproducible_enriched_re/", prefix, ".reproducible_enriched_re.tsv.gz")
)
