# Load required packages for data manipulation/plotting and label repulsion in ggplot2.
library(tidyverse)
library(ggrepel)

# Create output directories if they do not already exist (suppress warnings, create parents as needed).
dir.create("output/figures/secondary_figures/clip_scatter_re/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/enriched_re/", showWarnings = FALSE, recursive = TRUE)

# Parse command-line arguments: input tables and labels, plus output stem.
args = commandArgs(trailingOnly=TRUE)
re_data = read_tsv(args[1])  # Repeat-element count table per replicate (wide format with replicate columns).
re_annotations = read_tsv(args[2], col_names = c("chr","start","end","name","score","strand","repeat_name","repeat_class","repeat_family","gc"))  # Repeat annotations with GC content per locus.
model_data = read_tsv(args[3])  # Overdispersion model coefficients per GC bin or global fit.
input_replicate_label = args[4]  # Column name for the input replicate to evaluate.
clip_replicate_label = args[5]   # Column name for the CLIP replicate to evaluate.
output_stem = args[6]            # Basename used for plot and result outputs.

### Handle the case where re_data is empty by writing a dummy table and sentinel plot, then quitting.
if (length(re_data$repeat_name) < 1) {
    # Construct a dummy q_data-like table with the expected columns.
    dummy_q_data = tibble::tibble(
        chr            = NA_character_,
        start          = NA_real_,      # Use numeric NA for generality.
        end            = NA_real_,
        name           = NA_character_,
        score          = NA_real_,
        strand         = NA_character_,
        repeat_name    = NA_character_,
        repeat_class   = NA_character_,
        repeat_family  = NA_character_,
        gc             = NA_real_,
        gc_bin         = NA_character_,
        baseline_l2or  = NA_real_,
        clip           = NA_real_,
        input          = NA_real_,
        enrichment_l2or = NA_real_,
        pvalue         = NA_real_,
        simple         = NA_character_,
        qvalue         = NA_real_
    )

    # Write dummy enriched_re table for downstream steps.
    readr::write_tsv(
        dummy_q_data,
        paste0("output/secondary_results/enriched_re/", output_stem, ".enriched_re.tsv.gz")
    )

    # Create a sentinel plot indicating no data were available.
    pdf(
        paste0("output/figures/secondary_figures/clip_scatter_re/", output_stem, ".clip_test_distribution.pdf"),
        height = 1.8,
        width  = 2.8
    )
    par(mar = c(0,0,0,0))      # Critical fix: remove default margins.
    plot.new()
    text(0.5, 0.5, "No repeat-element data available for this sample.")
    dev.off()

    quit()
}

# Define helper link functions that use base-2 logs for convenience in l2or calculations.
logisticb2 = function(x) 1 / (1 + 2**-x)  # Inverse-logit on base-2 scale returning probability in (0,1).
logitb2 = function(x) log2(x / (1 - x))   # Logit on base-2 scale mapping probability to real line.

# Summarize GC content per repeat_name to later construct GC bins and baselines.
re_gc = re_annotations %>% group_by(repeat_name) %>% summarize(gc = mean(gc)) 

# Collapse overdispersion parameter to a single value (median on logit2 scale, then used later via inverse link).
model_overdispersion = median(logitb2(model_data$rho))

# Select and reshape the counts for the chosen input/CLIP replicates, remove other replicate columns, and compute GC-binned baselines.
selected_re_data = re_data %>% 
  rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)) %>% 
  select(-matches("(IP|IN)_[0-9]+$")) %>%                       # Drop any other replicate columns to avoid ambiguity.
  filter(input + clip > 0) %>%                                  # Keep loci with at least one read across input+CLIP.
  inner_join(re_gc) %>%                                         # Attach GC for each repeat_name.
  mutate(gc_bin = cut_number(gc,20)) %>%                        # Bin GC into 20 equal-count bins.
  group_by(gc_bin) %>%
  mutate(baseline_l2or = median( logitb2((clip + 1) / (clip + input + 2)) )) %>%  # Median log2-odds baseline per GC bin (with pseudocounts).
  ungroup

# Build per-locus test data: compute enrichment relative to GC-binned baseline and a beta-binomial p-value.
p_data = selected_re_data %>% 
  group_by(gc_bin, baseline_l2or, clip, input) %>%
  summarize %>% 
  mutate(
    enrichment_l2or = log2((clip + logisticb2(baseline_l2or)) / (input + 1 - logisticb2(baseline_l2or))) - baseline_l2or,  # L2OR enrichment above baseline.
    pvalue = pmax(1e-10, 1 - VGAM::pbetabinom(q = clip - 1, size = clip + input, prob = logisticb2(baseline_l2or), rho = model_overdispersion %>% logisticb2 ))  # Upper-tail beta-binomial p-value with floor.
  ) %>%
  inner_join(selected_re_data,.)  # Reattach metadata/annotation per locus.

# Multiple testing adjustment with an input+CLIP total-read threshold, and annotate a simple element category for plotting.
q_data = p_data %>% 
  mutate(simple = ifelse(grepl("\\)n",repeat_name), "kmer", "element")) %>%   # Tag k-mer rows by name pattern.
  group_by(above_threshold = input + clip >= 10) %>% 
  mutate(qvalue = ifelse(above_threshold, p.adjust(pvalue,"fdr"), NA)) %>%    # FDR only for sufficiently covered loci.
  ungroup %>% 
  arrange(pvalue, desc(enrichment_l2or))  # Order for consistent label selection downstream.

# Pick a few top labels (excluding k-mers) for annotation on the scatter plot.
label_data = q_data %>% filter(!grepl("\\)n$", repeat_name), enrichment_l2or > 2.5) %>% head(4)

# Symmetric y-limits based on the maximum absolute enrichment value.
enrichment_max = q_data$enrichment_l2or %>% abs %>% max

# Produce the scatter plot of enrichment vs total reads (log10 scale), coloring by significance and shaping by element class.
pdf(paste0('output/figures/secondary_figures/clip_scatter_re/', output_stem, '.clip_test_distribution.pdf'),height = 1.8, width = 2.8)
    q_data %>% replace_na(list(qvalue = 1)) %>% filter(clip + input >= 10) %>%
    ggplot(aes(clip + input, enrichment_l2or, color = ifelse(qvalue < 0.05, "q < 0.05", "NS"), shape = simple)) + 
        geom_hline(yintercept = 0) +
        geom_hline(linetype ="31", color ="#636363", yintercept = c(-2.5,2.5)) +
        theme_bw(base_size = 7) +
        scale_shape_manual(values = c(16,1)) +
        theme(legend.title=element_blank()) +
        expand_limits(y = c(-enrichment_max, enrichment_max)) +
        geom_point(size = 0.8, stroke = 0) + scale_x_log10() + xlab("Total reads") + ylab("Log2 enrichment") + 
        geom_text_repel(data = label_data, aes(label = repeat_name),min.segment.length = 0.1, color = "black",size = 2) + 
        theme(panel.grid.minor = element_blank()) + scale_color_manual(values = c("#bdbdbd","#f86808"))
        
dev.off()

# Write out the per-element statistics (with q-values and enrichment) for downstream use.
write_tsv(q_data %>% select(-above_threshold), paste0("output/secondary_results/enriched_re/", output_stem, ".enriched_re.tsv.gz"))