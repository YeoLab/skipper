# Load core tidyverse for data wrangling/IO/plotting and companion libs for t-SNE and label repulsion in scatter plots.
library(tidyverse)
library(Rtsne)
library(ggrepel)

# Create output directories if missing (idempotent) for tables and figures produced by this script.
dir.create("output/tsne/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/tsne/", showWarnings = FALSE, recursive = TRUE)

# Parse command-line arguments.
args = commandArgs(trailingOnly=TRUE)
win_directory = args[1]
ref_directory = args[2]
prefix = args[3]

# Discover all enriched window and repeat-element result files to be ingested.
enriched_window_files = list.files(path = win_directory, pattern = ".*\\.reproducible_enriched_windows.tsv.gz", full.names = TRUE)

# Read, keep nonempty, and bind enriched windows across experiments, tagging each row with its experiment label.
enriched_windows = enriched_window_files %>%
    setNames(sub("\\.reproducible_enriched_windows\\.tsv.gz", "", basename(.))) %>% 
    map(function(x) read_tsv(x)) %>% 
    Filter(function(x) nrow(x) > 0, .) %>%
    bind_rows(.id = "experiment_label") 

# Collect the union of all experiment IDs present in either input.
experiment_labels = enriched_windows$experiment_label %>% unique

# Load reference resources used to standardize features and to plot against ENCODE reference embeddings/clusters.
repeatmasker_labels = read_tsv(paste0(ref_directory, "/repeat_masker_hierarchy.tsv"))
reference_enrichments = read_tsv(paste0(ref_directory, "/encode3_eclip_enrichment.reference.tsv"))
reference_features = read_tsv(paste0(ref_directory, "/encode3_feature_summary.reference.tsv"))
reference_assignments = read_tsv(paste0(ref_directory, "/encode3_class_assignment.reference.tsv")) %>% select(-tsne_1,-tsne_2) %>%
	bind_rows(tibble(id = experiment_labels, rbp = NA, cluster = 0, class = "Query")) %>% 
	mutate(class = factor(class, levels = c("5' UTR","CDS","3' UTR","Splice site","Intron","MtRNA", "YRNA","snoRNA","tRNA/snRNA","Query")))

# Helpers to map detailed window/RE annotations to a simplified feature vocabulary for cross-sample comparison.
simplify_window_data = function(df) {
	df %>% mutate(
		feature = paste0(transcript_type_top,":",feature_type_top),
		feature = ifelse(grepl("Y_RNA|^RNY", gene_name), "Y_RNA", feature),
		feature = ifelse(grepl("^7SK$|^RN7SK", gene_name), "7SK", feature),
		feature = ifelse(chrom == "chrM", "Mt_RNA", feature),
		feature = ifelse(gene_type_top == "snRNA", "snRNA", feature),
		feature = ifelse(gene_type_top == "snoRNA", "snoRNA", feature),
		feature = ifelse(gene_type_top == "^rRNA", "rRNA", feature)
	)
}

# Count simplified feature occurrences per experiment for windows and REs, and merge with reference feature priors.
win_features = simplify_window_data(enriched_windows) %>% 
	group_by(id = experiment_label, feature) %>% count(name="value") %>% ungroup

# Join with reference background fractions, fill missing with priors, and aggregate per id/feature counts.
clip_count_data = win_features %>% group_by(id) %>% summarize(left_join(reference_features %>% mutate(id=unique(id)), tibble(.))) %>% 
	mutate(value = ifelse(is.na(value), global_fraction, value)) %>% group_by(id, feature, global_fraction) %>% summarize(clip_count = sum(value)) 

# Convert counts to fractions and compute information contribution (entropy_contribution) per feature and sample.
clip_enrichment_data = clip_count_data %>% group_by(id) %>% mutate(clip_fraction = clip_count / sum(clip_count), entropy_contribution = log2(clip_fraction / global_fraction) * clip_fraction) %>% arrange(entropy_contribution %>% desc)

# Build a correlation-distance matrix (1 − Pearson correlation) across all samples (reference + query) over the feature space.
enrichment_dist = bind_rows(reference_enrichments, clip_enrichment_data) %>% pivot_wider(values_from = "entropy_contribution",names_from="id",feature) %>% column_to_rownames("feature") %>% cor(m="p") %>% (function(x) 1 - x) %>% as.dist

# Run t-SNE on the distance matrix (is_distance=TRUE) for 2D embedding suitable for visualization.
set.seed(400)
tsne = Rtsne(enrichment_dist, perplexity = 20, is_distance = TRUE)

# Collect the embedding and join with reference class assignments, marking queries as “Query”.
tsne_data = tibble(id = enrichment_dist %>% as.matrix %>% row.names) %>% bind_cols(tsne$Y %>% as.data.frame %>% setNames(c("tsne_1","tsne_2"))) %>% inner_join(reference_assignments)

# Color palette for classes in the embedding plot.
samba_colors = c("#1B85ED", "#1AA2F6", "#00BFFF", "#4AC596", "#00CC00", "#A7D400", "#FFD700", "#FFBE00", "#FFA500", "#636363")

# Plot the t-SNE embedding, labeling query points, and save to PDF.
pdf(paste0("output/figures/tsne/", prefix, ".tsne_query.pdf"),height = 1.7, width = 2.6)
ggplot(tsne_data, aes(tsne_1, tsne_2,color = class)) + theme_bw(base_size = 7) + 
	geom_point(size = 0.8,stroke=0) + theme(panel.grid = element_blank(),aspect.ratio = 1,legend.key.size = unit(0.25,"cm"), legend.spacing.y = unit(0.2,"cm")) + 
	xlab("t-SNE 1") + ylab("t-SNE 2") + scale_color_manual(values = samba_colors) +
	geom_text_repel(data=tsne_data %>% filter(class == "Query"), aes(label = id), color="black", size = 1.4,min.segment.length = 0.1,segment.size = 0.2)
dev.off() 

# Persist the embedding coordinates and class labels for downstream analyses.
write_tsv(tsne_data, paste0("output/tsne/", prefix, ".tsne_query.tsv"))
