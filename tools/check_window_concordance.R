library(ggupset)
library(tidyverse)

dir.create("output/figures/enrichment_concordance/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/enrichment_reproducibility/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/enrichment_reproducibility/", showWarnings = FALSE, recursive = TRUE)
args = commandArgs(trailingOnly=TRUE)

data_directory = args[1]
prefix = args[2]

if(length(args) > 2) {
	blacklist = read_tsv(args[3], col_names = c("chr","start","end","name","score","strand"), col_types = "cddcdc")
} else {
	blacklist = tibble(chr=character(),start=numeric(),end=numeric(),name=character(),score=numeric(),strand=character())
}
tested_window_files = list.files(path = data_directory, pattern = paste0("^", prefix, "\\..*tested_windows.tsv.*"), full.names = TRUE)

tested_window_data = tested_window_files %>%
    setNames(sub("\\.tested_windows\\.tsv.*", "", basename(.))) %>% 
    map(function(x) read_tsv(x) %>% select(chr, start, end, name, score, strand, gc, gc_bin, qvalue) %>% mutate(name = as.character(name)) %>% anti_join(blacklist %>% select(-name))) %>% 
    bind_rows(.id = "clip_replicate_label") 

clip_replicate_labels = unique(tested_window_data$clip_replicate_label)
for(clip_replicate_1 in clip_replicate_labels) {
	for(clip_replicate_2 in clip_replicate_labels) {
		if(clip_replicate_1 < clip_replicate_2) {
			enriched_concordance_data = tested_window_data %>% 
				filter(clip_replicate_label %in% c(clip_replicate_1, clip_replicate_2)) %>% 
				pivot_wider(names_from = clip_replicate_label, values_from = qvalue) %>% 
				filter(!is.na(.data[[clip_replicate_1]]), !is.na(.data[[clip_replicate_2]])) %>% 
				group_by(replicate_1 = ifelse(.data[[clip_replicate_1]] < 0.2, "Enriched", "Not enriched"), replicate_2 = ifelse(.data[[clip_replicate_2]] < 0.2, "Enriched", "Not enriched")) %>% 
				count(name="count") %>% 
				mutate(replicate_1 = factor(replicate_1, levels = c("Not enriched", "Enriched"))) %>%
				mutate(replicate_2 = factor(replicate_2, levels = c("Enriched", "Not enriched")))
			
			if(nrow(enriched_concordance_data) < 4) {
				pdf(paste0("output/figures/enrichment_concordance/", prefix, ".enrichment_concordance.", clip_replicate_1, "_", clip_replicate_2, ".pdf"), height = 1, width = 2)
					print(ggplot() + annotate("text", x = 1, y = 1, label = "Insufficient data") + theme_void())
				dev.off()
				next	
			}
			
			odds_data = enriched_concordance_data %>% pivot_wider(names_from="replicate_2", values_from = "count") %>% 
				ungroup %>% select(-replicate_1) %>% fisher.test %>% (broom::tidy)	

			# Save odds_data to a TSV file
			output_tsv_path <- paste0("output/enrichment_reproducibility/", prefix, ".odds_data.tsv")
			write_tsv(odds_data, output_tsv_path)

			or_label = with(odds_data, ifelse(p.value < 0.05 & estimate > 1, paste0(sprintf(fmt="%.3g", odds_data$estimate),"x"), "NS"))
			enriched_mosaic_data = enriched_concordance_data %>% group_by(replicate_1) %>% mutate(width = sum(count), fraction = count/sum(count)) %>% ungroup %>% mutate(width = width / sum(width)*3.9)
			bar_break_y = enriched_mosaic_data %>% filter(replicate_1 == "Enriched", replicate_2 == "Not enriched") %>% pull(fraction)
			annotation_y = bar_break_y + 0.5*(1 - bar_break_y) 
			pdf(paste0("output/figures/enrichment_concordance/", prefix, ".enrichment_concordance.", clip_replicate_1, "_", clip_replicate_2, ".pdf"), height = 1.8, width = 3)
			print(
				enriched_mosaic_data %>%
				ggplot(aes(replicate_1, fraction * 1.95, fill = replicate_2, width = width)) + theme_minimal(base_size = 7) +
					geom_bar(stat = "identity",color = "black",size = 0.4, position = position_stack()) + theme(panel.grid = element_blank()) +  
					scale_fill_manual("", values = c("#8fb576", "#d9d9d9")) + coord_equal() +
					xlab(clip_replicate_1) + ylab(clip_replicate_2) + 
					scale_y_continuous(breaks = 1.95*c(0,0.25,0.5,0.75,1),labels = as.character(c(0,0.25,0.5,0.75,1))) +
					ggtitle(paste0(or_label, " odds ratio")) +
					# annotate("rect", xmin= "Enriched", xmax = y = 0.9, fill = "black", size = 2.8) +
					annotate("text", x= "Enriched", y = annotation_y * 1.95, label = "*", color = "white", size = 3.5)
			)
			dev.off()
		}
	}
}

# p_enriched = tested_window_data %>% group_by(S=qvalue < 0.2) %>% count %>% ungroup %>% with(n[S] / sum(n))

thresholded_window_data = tested_window_data %>% 
	group_by(chr, start, end, name, score, strand, status=ifelse(qvalue < 0.2, "Enriched", "Not enriched")) %>% 
	count(name="# Replicates")

enrichment_reproducibility_data = thresholded_window_data %>% 
	group_by(status, `# Replicates`) %>% count(name="N") %>% 
	pivot_wider(names_from=status,values_from = N,values_fill=0) %>%
	pivot_longer(names_to = "Status", values_to = "# Tested windows", - `# Replicates`)

pdf(paste0("output/figures/enrichment_reproducibility/", prefix, ".enrichment_reproducibility.pdf"), height = 1.2, width = 2)
ggplot(enrichment_reproducibility_data, aes(`# Replicates`, `# Tested windows`, color = Status)) + theme_bw(base_size=7)+
	geom_point() + scale_y_log10(limits = c(1, NA)) + scale_color_manual("", values = c("#fec332", "#bdbdbd")) +
	theme(panel.grid.minor = element_blank()) + scale_x_continuous(breaks = seq(1, length(tested_window_files)))
dev.off()

write_tsv(enrichment_reproducibility_data, paste0("output/enrichment_reproducibility/", prefix, ".enrichment_reproducibility.tsv"))
