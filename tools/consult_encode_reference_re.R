# Load core tidyverse for data wrangling/IO/plotting and companion libs for t-SNE and label repulsion in scatter plots.
library(tidyverse)
library(Rtsne)
library(ggrepel)

# Create output directories if missing (idempotent) for tables and figures produced by this script.
dir.create("output/tsne_re/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/tsne_re/", showWarnings = FALSE, recursive = TRUE)

# Parse command-line arguments.
args = commandArgs(trailingOnly=TRUE)
re_directory = args[1]
ref_directory = args[2]
prefix = args[3]

# Discover all enriched window and repeat-element result files to be ingested.
enriched_re_files = list.files(path = re_directory, pattern = ".*\\.reproducible_enriched_re.tsv.gz", full.names = TRUE)

# Same for enriched repeat elements (REs).
enriched_re = enriched_re_files %>%
    setNames(sub("\\.reproducible_enriched_re\\.tsv.gz", "", basename(.))) %>% 
    map(function(x) read_tsv(x)) %>% 
    Filter(function(x) nrow(x) > 0, .) %>%
    bind_rows(.id = "experiment_label") 

# Collect the union of all experiment IDs present in either input.
experiment_labels = enriched_re$experiment_label %>% unique

# Load reference resources used to standardize features and to plot against ENCODE reference embeddings/clusters.
repeatmasker_labels = read_tsv(paste0(ref_directory, "/repeat_masker_hierarchy.tsv"))
reference_enrichments = read_tsv(paste0(ref_directory, "/encode3_eclip_enrichment.reference.tsv"))
reference_features = read_tsv(paste0(ref_directory, "/encode3_feature_summary.reference.tsv"))
reference_assignments = read_tsv(paste0(ref_directory, "/encode3_class_assignment.reference.tsv")) %>% select(-tsne_1,-tsne_2) %>%
	bind_rows(tibble(id = experiment_labels, rbp = NA, cluster = 0, class = "Query")) %>% 
	mutate(class = factor(class, levels = c("5' UTR","CDS","3' UTR","Splice site","Intron","MtRNA", "YRNA","snoRNA","tRNA/snRNA","Query")))

# Helpers to map detailed window/RE annotations to a simplified feature vocabulary for cross-sample comparison.
simplify_re_data = function(df) {
	filter(df, simple == "element", enrichment_l2or_mean > 2.5) %>% mutate(
		feature = replace(repeat_name, repeat_class == "snRNA", "snRNA"),
		feature = replace(feature, grepl("^HY[0-9]$", repeat_name), "Y_RNA"),
		feature = replace(feature, repeat_class == "tRNA", "tRNA"),
		feature = replace(feature, repeat_family == "L1", "L1"),
		feature = replace(feature, repeat_family == "Alu", "Alu"),
		feature = replace(feature, repeat_family == "L1_AS", "Antisense_L1"),
		feature = replace(feature, repeat_family == "Alu_AS", "Antisense_Alu"),
		feature = replace(feature, repeat_class == "LTR", "LTR"),
		feature = replace(feature, repeat_name == "7SK", "7SK"),
		feature = replace(feature, repeat_class == "LINE" & repeat_family != "L1", "Misc_LINE"),
		feature = replace(feature, !feature %in% c("snRNA","Y_RNA","tRNA","L1","Alu","Antisense_L1","Antisense_Alu","LTR","7SK","Misc_LINE"), "Other_RE")
	)
}

# Count simplified feature occurrences per experiment for windows and REs, and merge with reference feature priors.
re_features = inner_join(enriched_re,repeatmasker_labels) %>% simplify_re_data %>% 
	group_by(id = experiment_label, feature) %>% count(name = "value")  %>% ungroup

# Join with reference background fractions, fill missing with priors, and aggregate per id/feature counts.
clip_count_data = re_features %>% group_by(id) %>% summarize(left_join(reference_features %>% mutate(id=unique(id)), tibble(.))) %>% 
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
pdf(paste0("output/figures/tsne_re/", prefix, ".tsne_re_query.pdf"),height = 1.7, width = 2.6)
ggplot(tsne_data, aes(tsne_1, tsne_2,color = class)) + theme_bw(base_size = 7) + 
	geom_point(size = 0.8,stroke=0) + theme(panel.grid = element_blank(),aspect.ratio = 1,legend.key.size = unit(0.25,"cm"), legend.spacing.y = unit(0.2,"cm")) + 
	xlab("t-SNE 1") + ylab("t-SNE 2") + scale_color_manual(values = samba_colors) +
	geom_text_repel(data=tsne_data %>% filter(class == "Query"), aes(label = id), color="black", size = 1.4,min.segment.length = 0.1,segment.size = 0.2)
dev.off() 

# Persist the embedding coordinates and class labels for downstream analyses.
write_tsv(tsne_data, paste0("output/tsne_re/", prefix, ".tsne_re_query.tsv"))
