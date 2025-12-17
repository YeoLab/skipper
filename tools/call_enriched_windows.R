# Load plotting palettes and tidyverse for data manipulation/plotting, readr, dplyr, ggplot2, etc.
library(tidyverse)
library(viridis)

# Parse CLI args and create all output directories (idempotent). 
args = commandArgs(trailingOnly=TRUE)
dir.create("output/secondary_results/threshold_scan/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/tested_windows/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/enriched_windows/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/enrichment_summaries/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/all_reads/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/secondary_figures/threshold_scan/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/secondary_figures/enriched_windows/", showWarnings = FALSE)
dir.create("output/figures/secondary_figures/all_reads/", showWarnings = FALSE)

# Read inputs: per-window counts, accession ranking table, feature annotations, and model parameters, plus labels/stem. 
count_data = read_tsv(args[1])
accession_data = read_tsv(args[2]) %>% arrange(rank)
feature_annotations = read_tsv(args[3])
model_data = read_tsv(args[4])
input_replicate_label = args[5]
clip_replicate_label = args[6]
output_stem = args[7]

# Convenience link functions on base-2 scale for logits and inverse-logits. 
logisticb2 = function(x) 1 / (1 + 2**-x)
logitb2 = function(x) log2(x / (1 - x))

# Prepare exon subtype ordering and plotting orders for features and transcript types. 
exon_subtypes = accession_data$exon_subtype %>% unique
protein_coding_subtype = accession_data$exon_subtype[accession_data$accession == "protein_coding"] %>% head(1)
prioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) < 1]
unprioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) >= 1]
feature_plot_order = c("CDS_SOLITARY", "CDS_START","CDS_STOP","CDS","UTR5","UTR3",paste0("EXON_", prioritized_exon_subtypes),paste0("EXON_", unprioritized_exon_subtypes),"SSB_ADJ","SSB_PROX","SS3_ADJ","SS3_PROX","SS5_ADJ","SS5_PROX","PRIMIRNA","INTRON") %>% rev
transcript_plot_order = feature_annotations %>% group_by(transcript_type_top) %>% count(sort=TRUE) %>% mutate(tname = gsub("_","\n",transcript_type_top)) %>% pull(tname) 

# Keep windows with any reads and form GC deciles within which to compute baselines. 
count_gc_data = count_data[select(count_data, matches("(IP|IN)_[0-9]+$")) %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(gc,10)) %>% filter(.data[[clip_replicate_label]] + .data[[input_replicate_label]] > 0)

# Rename chosen replicates to `input`/`clip`, compute GC-bin baselines (mean clip fraction), and attach back to rows. 
processed_count_data = count_gc_data %>% rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)) %>% 
	summarize(baseline_l2or = mean((clip / (clip + input))) %>% logitb2) %>% 
	inner_join(select(count_gc_data %>% rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)), -matches("(IP|IN)_[0-9]+$")),.)

# Set beta-binomial overdispersion (rho) from model, on probability scale. 
model_overdispersion = VGAM::logitlink(median(model_data$rho), inverse=TRUE)

# Collapse to unique (clip,input,gc_bin,baseline) combos, compute enrichment and p-values, and rejoin to per-window rows. 
p_data = processed_count_data %>% group_by(clip, input, gc_bin, baseline_l2or) %>%
	summarize %>% mutate(enrichment_l2or = log2((clip + logisticb2(baseline_l2or)) / (input + 1 - logisticb2(baseline_l2or))) - baseline_l2or) %>%
	mutate(pvalue = pmax(1e-12, 1 - VGAM::pbetabinom(q = clip - 1, size = clip + input, prob = logisticb2(baseline_l2or), rho = model_overdispersion))) %>%
	inner_join(processed_count_data,.)

# Convenience link functions on base-2 scale for logits and inverse-logits. 
logisticb2 = function(x) 1 / (1 + 2**-x)
logitb2 = function(x) log2(x / (1 - x))

# Prepare exon subtype ordering and plotting orders for features and transcript types. 
exon_subtypes = accession_data$exon_subtype %>% unique
protein_coding_subtype = accession_data$exon_subtype[accession_data$accession == "protein_coding"] %>% head(1)
prioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) < 1]
unprioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) >= 1]
feature_plot_order = c("CDS_SOLITARY", "CDS_START","CDS_STOP","CDS","UTR5","UTR3",paste0("EXON_", prioritized_exon_subtypes),paste0("EXON_", unprioritized_exon_subtypes),"SSB_ADJ","SSB_PROX","SS3_ADJ","SS3_PROX","SS5_ADJ","SS5_PROX","PRIMIRNA","INTRON") %>% rev
transcript_plot_order = feature_annotations %>% group_by(transcript_type_top) %>% count(sort=TRUE) %>% mutate(tname = gsub("_","\n",transcript_type_top)) %>% pull(tname) 

# Keep windows with any reads and form GC deciles within which to compute baselines. 
count_gc_data = count_data[select(count_data, matches("(IP|IN)_[0-9]+$")) %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(gc,10)) %>% filter(.data[[clip_replicate_label]] + .data[[input_replicate_label]] > 0)

# Rename chosen replicates to `input`/`clip`, compute GC-bin baselines (mean clip fraction), and attach back to rows. 
processed_count_data = count_gc_data %>% rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)) %>% 
	summarize(baseline_l2or = mean((clip / (clip + input))) %>% logitb2) %>% 
	inner_join(select(count_gc_data %>% rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)), -matches("(IP|IN)_[0-9]+$")),.)

# Set beta-binomial overdispersion (rho) from model, on probability scale. 
model_overdispersion = VGAM::logitlink(median(model_data$rho), inverse=TRUE)

# Collapse to unique (clip,input,gc_bin,baseline) combos, compute enrichment and p-values, and rejoin to per-window rows. 
p_data = processed_count_data %>% group_by(clip, input, gc_bin, baseline_l2or) %>%
	summarize %>% mutate(enrichment_l2or = log2((clip + logisticb2(baseline_l2or)) / (input + 1 - logisticb2(baseline_l2or))) - baseline_l2or) %>%
	mutate(pvalue = pmax(1e-12, 1 - VGAM::pbetabinom(q = clip - 1, size = clip + input, prob = logisticb2(baseline_l2or), rho = model_overdispersion))) %>%
	inner_join(processed_count_data,.)

# Global IP fraction across all windows, used in some heuristics below. 
p_clip = with(p_data, sum(clip) / sum(clip + input))

# Add total counts (if not already present).
p_data <- p_data %>%
  mutate(total_counts = input + clip)

# Figure out an upper bound for the threshold scan.
threshold_max = p_data %>% mutate(total_counts = input + clip) %>% arrange(desc(total_counts)) %>% head(100) %>% tail(1) %>% pull(total_counts)

# Cap at 500 and ensure it's at least 2.
threshold_max <- min(500, threshold_max)

thresholds <- 2:threshold_max

threshold_data <- tibble(
threshold  = thresholds,
n_enriched = vapply(
    thresholds,
    function(th) {
        idx <- p_data$total_counts >= th
        if (!any(idx)) return(0L)
        padj_th <- p.adjust(p_data$pvalue[idx], method = "fdr")
        sum(padj_th < 0.2)
    },
    integer(1)
)
)

write_tsv(
  threshold_data,
  paste0("output/secondary_results/threshold_scan/", output_stem, ".threshold_data.tsv")
)

# Choose the threshold with the maximum number of discoveries. 
optimized_threshold = threshold_data %>% arrange(n_enriched %>% desc) %>% head(1) %>% pull(threshold)

p_data$total_counts = NULL

# Compute q-values using the optimized threshold and carry an indicator for filtering/outputs. 
q_data = p_data %>% group_by(above_threshold = input + clip >= optimized_threshold) %>% 
	mutate(qvalue = ifelse(above_threshold, p.adjust(pvalue,"fdr"), NA)) %>% ungroup

# Summarize total read fractions by feature class for input vs clip and plot bars. 
all_reads_fractions_feature_data = q_data %>% inner_join(select(feature_annotations, name, feature_type_top)) %>% 
	group_by(feature_type_top) %>% summarize(input = sum(input), clip = sum(clip)) %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>%
	mutate(input = input / sum(input), clip = clip / sum(clip)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev))

pdf(paste0('output/figures/secondary_figures/all_reads/', output_stem, '.all_reads_fractions.feature.pdf'), height = 1.8, width = 2.8)
all_reads_fractions_feature_data %>%
	pivot_longer(names_to = "replicate", values_to = "fraction",-c(feature_type_top, feature_group)) %>%
	mutate(replicate = factor(replicate, levels = c("input", "clip"))) %>%
	mutate(clip_replicate = clip_replicate_label) %>%
ggplot(aes(x=feature_type_top, fraction, fill = feature_group, group = replicate)) + theme_bw(base_size = 7) + 
	geom_bar(stat = "identity",position = "dodge",size = 0.2) + 
	geom_bar(data = mutate(all_reads_fractions_feature_data %>% pivot_longer(names_to = "replicate", values_to = "fraction",-c(feature_type_top, feature_group)), replicate = factor(replicate, levels = c("input", "clip")), fraction = ifelse(replicate == "clip", 0, fraction)), stat = "identity",position = "dodge",size = 0.2, fill = "#ffffff", alpha = 0.5) + 
	scale_color_manual(values = c("white","black")) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
	ylab("Fraction of reads") + xlab("Type of feature") + 
	facet_wrap(~clip_replicate)
dev.off()

# Count all windows per feature type for denominators in rate/odds summaries. 
all_window_feature_data = feature_annotations %>% group_by(feature_type_top) %>% count(name = "n_windows") %>% ungroup

# Persist the tested-windows table (above threshold) with formatted numeric columns. 
tested_window_data = q_data %>% filter(above_threshold) %>% select(-above_threshold) %>% mutate(across(c("baseline_l2or", "enrichment_l2or","pvalue","qvalue"), ~ sprintf(.x, fmt = "%.6g")))
write_tsv(tested_window_data, paste0("output/secondary_results/tested_windows/", output_stem, ".tested_windows.tsv.gz"))

# Extract enriched windows (q < 0.2), join annotations, and sort by q-value. 
enriched_window_data = q_data %>% filter(above_threshold, qvalue < 0.2) %>% select(-above_threshold) %>% left_join(feature_annotations) %>% arrange(qvalue)
write_tsv(enriched_window_data %>% mutate(across(c("baseline_l2or", "enrichment_l2or","pvalue","qvalue"), ~ sprintf(.x, fmt = "%.6g"))), paste0("output/secondary_results/enriched_windows/", output_stem, ".enriched_windows.tsv.gz"))

# If no enriched windows, emit placeholder PDFs and empty summary TSVs, then quit gracefully. 
if(nrow(enriched_window_data) == 0) {
	for(output_file in c(
		paste0('output/figures/secondary_figures/enriched_windows/', output_stem, '.enriched_window_rates.pdf'),
		paste0('output/figures/secondary_figures/enriched_windows/', output_stem, '.enriched_window_counts.linear.pdf'),
		paste0('output/figures/secondary_figures/enriched_windows/', output_stem, '.enriched_window_counts.per_gene_feature.pdf')
	))
	{
		pdf(output_file, height = 1, width = 2)
			print(ggplot() + annotate("text", x = 1, y = 1, label = "No enriched windows") + theme_void())
		dev.off()
	}
	file.create(paste0("output/secondary_results/enrichment_summaries/", output_stem, ".enriched_window_feature_summary.tsv"))
	file.create(paste0("output/secondary_results/enrichment_summaries/", output_stem, ".enriched_window_transcript_summary.tsv"))
	file.create(paste0("output/secondary_results/enrichment_summaries/", output_stem, ".enriched_window_gene_summary.tsv"))
	quit()
}

# Summarize window testing and enrichment rates per feature type. 
enriched_window_feature_summary = left_join(
		all_window_feature_data,
		tested_window_data %>% inner_join(select(feature_annotations, name, feature_type_top)) %>% group_by(feature_type_top) %>% count(name = "n_tested")) %>%
	left_join(
		enriched_window_data %>% group_by(feature_type_top) %>% count(name = "n_enriched")
	) %>% replace_na(list(n_tested = 0, n_enriched = 0)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	mutate(`Test rate` = n_tested / n_windows, `Positive rate` = n_enriched / n_tested, `Enrichment rate` = n_enriched / n_windows)

# Plot enriched window counts per feature (linear scale). 
pdf(paste0('output/figures/secondary_figures/enriched_windows/', output_stem, '.enriched_window_counts.linear.pdf'), height = 1.8, width = 2.2)
enriched_window_feature_summary %>% mutate(replicate = clip_replicate_label) %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, n_enriched, fill = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_bar(stat="identity") + facet_grid(.~replicate) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched")
dev.off()

# Save feature-level enrichment summary table. 
write_tsv(enriched_window_feature_summary, paste0("output/secondary_results/enrichment_summaries/", output_stem, ".enriched_window_feature_summary.tsv"))

# Build transcript-by-feature window counts and tested/enriched summaries. 
all_window_transcript_data = feature_annotations %>% group_by(transcript_type_top, feature_type_top) %>% count(name = "n_windows") %>% ungroup

enriched_window_transcript_summary = left_join(
		all_window_transcript_data,
		tested_window_data %>% inner_join(select(feature_annotations, name, transcript_type_top, feature_type_top)) %>% group_by(transcript_type_top, feature_type_top) %>% count(name = "n_tested")) %>%
		replace_na(list(n_windows = 0, n_tested = 0)) %>%
	left_join(
		enriched_window_data %>% group_by(transcript_type_top, feature_type_top) %>% count(name = "n_enriched")
	) %>% replace_na(list(n_enriched = 0)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	mutate(transcript_type_top = factor(gsub("_", "\n", transcript_type_top), levels = transcript_plot_order %>% rev)) %>%
	mutate(`Test rate` = n_tested / n_windows, `Positive rate` = n_enriched / n_tested, `Enrichment rate` = n_enriched / n_windows)

# Pseudocount-based odds ratio across transcript-feature combinations. 
enriched_window_transcript_summary$pseudocount = enriched_window_transcript_summary %>% with(nrow(.) * (n_windows / (sum(n_windows))))
enriched_window_transcript_summary$odds_ratio = enriched_window_transcript_summary %>% with((n_enriched + pseudocount) / (sum((n_enriched + pseudocount)) - (n_enriched + pseudocount)) / (n_windows / (sum(n_windows) - n_windows)))

# Save transcript-feature enrichment summary. 
write_tsv(enriched_window_transcript_summary %>% mutate(transcript_type_top = gsub("\n", "_", transcript_type_top)), paste0("output/secondary_results/enrichment_summaries/", output_stem, ".enriched_window_transcript_summary.tsv"))

# Summarize enriched windows by gene and feature, bucket counts, and plot stacked bars. 
enriched_window_gene_summary = enriched_window_data %>% group_by(gene_name,feature_type_top) %>% count(name = "n_enriched") %>% 
		group_by(feature_type_top, n_enriched) %>% mutate(n_enriched = pmin(n_enriched, 5)) %>% count(name="n_genes") %>% 
		mutate(n_enriched = ifelse(n_enriched <5, as.character(n_enriched), "5+"))

pdf(paste0('output/figures/secondary_figures/enriched_windows/', output_stem, '.enriched_window_counts.per_gene_feature.pdf'), height = 2.2, width = 2.8)
enriched_window_gene_summary %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	mutate(replicate = clip_replicate_label) %>%
ggplot(aes(feature_type_top, n_genes, fill = n_enriched)) + theme_bw(base_size = 7) + 
	geom_bar(stat = "identity", position="stack") + facet_wrap(~replicate, nrow = 1) + 
	scale_fill_viridis("# enriched\nwindows",discrete=TRUE, option="plasma") +
	ylab("# genes") + xlab("Type of feature") + 
	theme(legend.position = "top", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank())
dev.off()

# Save gene-feature enrichment summary table. 
write_tsv(enriched_window_gene_summary, paste0("output/secondary_results/enrichment_summaries/", output_stem, ".enriched_window_gene_summary.tsv"))

