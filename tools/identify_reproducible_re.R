library(tidyverse)

dir.create("output/reproducible_enriched_re/", showWarnings = FALSE, recursive = TRUE)
args = commandArgs(trailingOnly=TRUE)

data_directory = args[1]
prefix = args[2]

enriched_re_files = list.files(path = data_directory, pattern = paste0("^", prefix, ".*enriched_re.tsv.gz"), full.names = TRUE)
enriched_re_data = enriched_re_files %>%
    setNames(sub("\\.enriched_windows\\.tsv.gz", "", basename(.))) %>% 
    map(function(x) read_tsv(x)) %>% 
    Filter(function(x) nrow(x) > 0, .) %>%
    bind_rows(.id = "clip_replicate_label") 

reproducible_enriched_re_data = enriched_re_data %>% group_by(across(-c(clip_replicate_label,gc_bin,baseline_l2or,input,clip,enrichment_l2or,pvalue,qvalue))) %>%
	summarize(
		input_sum = sum(input),
		clip_sum = sum(clip),
		enrichment_n = sum(qvalue < 0.05),
		enrichment_l2or_min = min(enrichment_l2or),
		enrichment_l2or_mean = mean(enrichment_l2or),
		enrichment_l2or_max = max(enrichment_l2or),
		p_max = max(pvalue),
		p_min = min(pvalue),
		q_max = max(qvalue),
		q_min = min(qvalue)
	) %>% filter(enrichment_n > 1) %>% arrange(p_min,enrichment_l2or_mean %>% desc)
	
write_tsv(reproducible_enriched_re_data, paste0("output/reproducible_enriched_re/", prefix, ".reproducible_enriched_re.tsv.gz"))
