library(tidyverse)


# Command line arguments.
args = commandArgs(trailingOnly=TRUE)
enriched_window_file = args[1]   # TSV of enriched windows across experiments
term_gmt_file = args[2]          # Gene set terms (GMT format)
prefix = args[3]                 # Output file prefix (e.g. "encode3_go_terms")


# Load input data.
encode_windows = read_tsv(enriched_window_file)
terms = fgsea::gmtPathways(term_gmt_file)

# Define gene universe: all genes present in both GMT terms and enriched windows.
gene_universe = terms %>% 
  unlist %>% unique %>% 
  intersect(encode_windows$gene_name %>% strsplit(split=":") %>% unlist %>% unique)

# Keep only valid terms:
# - at least 10 genes
# - fewer than 200 genes
valid_terms = lapply(terms, function(term_list) intersect(term_list, gene_universe)) %>%
  Filter(f=function(x) (length(x) >= 10) & (length(x) < 200))


# Aggregate enriched window counts per gene and per experiment.
encode_clip_gene_data = encode_windows %>%
  group_by(id) %>%
  summarize(table(gene = strsplit(gene_name,split=":") %>% unlist) %>% as_tibble) %>%
  filter(gene %in% gene_universe) %>%
  rename(n_windows_enriched = n) %>%
  group_by(id) %>%
  mutate(n_total_windows_clip = sum(n_windows_enriched))

# Total enriched windows across all experiments.
total_enriched_windows = encode_clip_gene_data %>% with(sum(n_windows_enriched))


# Calculate term sizes and enrichment counts.
# For each term, count how many genes and windows are enriched.
encode_term_sizes = lapply(valid_terms, function(term_list) length(term_list)) %>%
  as_tibble %>%
  pivot_longer(all_of(names(.)),names_to="term",values_to="n_genes_term") 

encode_term_data = lapply(valid_terms, function(term_list) {
    filter(encode_clip_gene_data, gene %in% term_list) %>%
      group_by(id) %>%
      summarize(n_genes_enriched = n(), n_windows_enriched = sum(n_windows_enriched))
  }) %>%
  bind_rows(.id = "term") %>%
  group_by(id) %>%
  summarize(left_join(encode_term_sizes %>% mutate(id = unique(id)),.)) %>%
  mutate(
    n_genes_enriched = ifelse(is.na(n_genes_enriched), 0, n_genes_enriched),
    n_windows_enriched = ifelse(is.na(n_windows_enriched), 0, n_windows_enriched)
  ) %>%
  group_by(term) %>%
  mutate(n_total_windows_term = sum(n_windows_enriched)) %>%
  inner_join(encode_clip_gene_data %>% distinct(id,n_total_windows_clip)) %>%
  ungroup


# Save reference file: background expected windows per term.
encode_term_data %>%
  distinct(term,n_term_windows=n_total_windows_term, f_term_windows = n_total_windows_term / total_enriched_windows) %>%
  write_tsv(paste0(prefix, ".reference.tsv.gz"))


# Binomial test per term/experiment.
# Construct contingency table counts:
# a = total background windows outside both sets
# b = total windows in experiment but not enriched in term
# c = total windows in term but not enriched in experiment
# d = overlap (windows enriched in both experiment and term)
# Then test enrichment with binomial test
binom_data = encode_term_data %>%
  mutate(
    a = total_enriched_windows - n_total_windows_clip - n_total_windows_clip + n_windows_enriched,
    b = n_total_windows_clip - n_windows_enriched,
    c = n_total_windows_term - n_windows_enriched,
    d = n_windows_enriched,
    p = c/a,         # background probability
    x = d,           # successes (observed overlap)
    n = b+d          # total trials (all experiment windows)
  ) %>%
  rowwise %>%
  summarize(
    (broom::tidy)(binom.test(x= x, n=n,p=p)) %>% 
      ungroup %>% 
      mutate(id=id,term=term)
  ) %>%
  inner_join(encode_term_data)


# Save binomial test results.
write_tsv(binom_data,paste0(prefix, ".results.tsv.gz"))