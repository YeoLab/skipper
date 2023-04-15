library(tidyverse)

dir.create("output/figures/reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)
args = commandArgs(trailingOnly=TRUE)

data_directory = args[1]
prefix = args[2]
if(length(args) > 2) {
	blacklist = read_tsv(args[3], col_names = c("chr","start","end","name","score","strand"), col_types = "cddcdc")
} else {
	blacklist = tibble(chr=character(),start=numeric(),end=numeric(),name=character(),score=numeric(),strand=character())
}
enriched_window_files = list.files(path = data_directory, pattern = paste0("^", prefix, ".*enriched_windows.tsv.gz"), full.names = TRUE)

# specify col_type to handle case of no enriched windows
enriched_window_data = enriched_window_files %>%
    setNames(sub("\\.enriched_windows\\.tsv.gz", "", basename(.))) %>% 
    map(function(x) read_tsv(x, col_types="cddcdcdddcddddcddccccccccc") %>% mutate(name = as.character(name)) %>% anti_join(blacklist %>% select(-name))) %>% 
    Filter(function(x) nrow(x) > 0, .) %>%
    bind_rows(.id = "clip_replicate_label") 

if(nrow(enriched_window_data %>% group_by(name) %>% filter(n() > 1)) == 0) {
	pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.linear.pdf'), height = 1, width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "No data") + theme_void())
	dev.off()
	pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.log10.pdf'), height = 1, width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "No data") + theme_void())
	dev.off()
	reproducible_enriched_window_data = enriched_window_data %>% group_by(across(-c(clip_replicate_label,baseline_l2or,input,clip,enrichment_l2or,pvalue,qvalue))) %>% summarize %>% head(0)
	write_tsv(reproducible_enriched_window_data, paste0("output/reproducible_enriched_windows/", prefix, ".reproducible_enriched_windows.tsv.gz"))
	quit()
}	
reproducible_enriched_window_data = enriched_window_data %>% group_by(across(-c(clip_replicate_label,baseline_l2or,input,clip,enrichment_l2or,pvalue,qvalue))) %>%
	summarize(
		input_sum = sum(input),
		clip_sum = sum(clip),
		enrichment_n = sum(qvalue < 0.2),
		enrichment_l2or_min = min(enrichment_l2or),
		enrichment_l2or_mean = mean(enrichment_l2or),
		enrichment_l2or_max = max(enrichment_l2or),
		p_max = max(pvalue),
		p_min = min(pvalue),
		q_max = max(qvalue),
		q_min = min(qvalue)
	) %>% filter(enrichment_n > 1) %>% arrange(enrichment_l2or_mean %>% desc)

pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.linear.pdf'), height = 1.8, width = 2.2)
reproducible_enriched_window_data %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, fill = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_bar() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched")
dev.off()

pdf(paste0('output/figures/reproducible_enriched_windows/', prefix, '.reproducible_enriched_window_counts.log10.pdf'), height = 1.8, width = 2.2)
reproducible_enriched_window_data %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% group_by(feature_group, feature_type_top) %>% count %>%
ggplot(aes(feature_type_top, n, color = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_point(stroke = 0) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched") + scale_y_log10()
dev.off()

write_tsv(reproducible_enriched_window_data, paste0("output/reproducible_enriched_windows/", prefix, ".reproducible_enriched_windows.tsv.gz"))
