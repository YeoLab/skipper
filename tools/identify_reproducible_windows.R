library(tidyverse)

# Create output directories for figures and data.
dir.create("output/figures/secondary_figures/unfiltered_reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/unfiltered_reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)

# Command-line arguments:
args = commandArgs(trailingOnly=TRUE)
data_directory = args[1]
prefix = args[2]
repro_def = args[3]

print(repro_def)

# Collect all enriched window files for the given experiment prefix.
enriched_window_files = list.files(
  path = data_directory,
  pattern = paste0("^", prefix, "\\..*enriched_windows.tsv.gz"),
  full.names = TRUE
)

if (repro_def == "ALL") {
    cutoff = length(enriched_window_files)
} else {
    cutoff = as.numeric(repro_def)
}

enriched_window_schema = readr::read_tsv(
  enriched_window_files[[1]],
  col_types = "cddcdcdddcddddcddccccccccc"
) %>%
  slice(0) %>%
  mutate(clip_replicate_label = character())

enriched_window_data = enriched_window_files %>%
  setNames(sub("\\.enriched_windows\\.tsv.gz", "", basename(.))) %>%
  map(\(x) {
    read_tsv(x, col_types = "cddcdcdddcddddcddccccccccc") %>%
      mutate(name = as.character(name))
  }) %>%
  purrr::keep(\(x) nrow(x) > 0) %>%
  bind_rows(.id = "clip_replicate_label") %>%
  { bind_rows(enriched_window_schema, .) }

# Force numeric types.
enriched_window_data = enriched_window_data %>%
  mutate(across(
    c(input, clip, enrichment_l2or, pvalue, qvalue),
    ~ suppressWarnings(as.numeric(.))
  ))

# Handle case: no enriched windows across all replicates.
if (nrow(enriched_window_data) == 0){
	# Construct an empty dataframe with expected columns
	columns= c("chr","start","end","name","score","strand","gc",
	"gc_bin","chrom","feature_id","feature_bin","feature_type_top","feature_types",
	"gene_name","gene_id","transcript_ids","gene_type_top","transcript_type_top",
	"gene_types","transcript_types", "input_sum","clip_sum","enrichment_n",
	"enrichment_l2or_min","enrichment_l2or_mean","enrichment_l2or_max","p_max","p_min",
	"q_max","q_min") 
	reproducible_enriched_window_data = data.frame(matrix(nrow = 0, ncol = length(columns))) 
	colnames(reproducible_enriched_window_data) = columns

	# Save empty table and exit
	write_tsv(reproducible_enriched_window_data, paste0("output/secondary_results/unfiltered_reproducible_enriched_windows/", prefix, ".unfiltered_reproducible_enriched_windows.tsv.gz"))
	quit()
}

# Handle case: only single-replicate enrichment (no overlap).
if (nrow(enriched_window_data %>% group_by(name) %>% filter(n() > 1)) == 0) {
	# Construct an empty dataframe with expected columns
	columns= c("chr","start","end","name","score","strand","gc",
	"gc_bin","chrom","feature_id","feature_bin","feature_type_top","feature_types",
	"gene_name","gene_id","transcript_ids","gene_type_top","transcript_type_top",
	"gene_types","transcript_types", "input_sum","clip_sum","enrichment_n",
	"enrichment_l2or_min","enrichment_l2or_mean","enrichment_l2or_max","p_max","p_min",
	"q_max","q_min") 
	reproducible_enriched_window_data = data.frame(matrix(nrow = 0, ncol = length(columns))) 
	colnames(reproducible_enriched_window_data) = columns

	# Save empty table and exit
	write_tsv(reproducible_enriched_window_data, paste0("output/secondary_results/unfiltered_reproducible_enriched_windows/", prefix, ".unfiltered_reproducible_enriched_windows.tsv.gz"))
	quit()
}	

# Prevents an annoying bug. 
stopifnot(is.numeric(enriched_window_data$qvalue))
           
# Aggregate reproducible enriched windows across replicates.
reproducible_enriched_window_data = enriched_window_data %>%
	group_by(chr,start,end,name,score,strand,gc,gc_bin,feature_id,feature_bin,chrom,feature_type_top,feature_types,
            gene_name, gene_id, transcript_ids, gene_type_top, transcript_type_top, gene_types, transcript_types) %>%
	summarize(input_sum = sum(input), clip_sum = sum(clip), enrichment_n = sum(qvalue < 0.2),
              enrichment_l2or_min = min(enrichment_l2or), enrichment_l2or_mean = mean(enrichment_l2or), 
              enrichment_l2or_max = max(enrichment_l2or), p_max = max(pvalue), p_min = min(pvalue),
              q_max = max(qvalue), q_min = min(qvalue)                                  
	) %>%
	filter(enrichment_n >= cutoff) %>%
	arrange(desc(enrichment_l2or_mean))

# Save reproducible enriched window data.
write_tsv(reproducible_enriched_window_data, paste0("output/secondary_results/unfiltered_reproducible_enriched_windows/", prefix, ".unfiltered_reproducible_enriched_windows.tsv.gz"))
