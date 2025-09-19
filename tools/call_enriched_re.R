library(tidyverse)
library(ggrepel)

dir.create("output/figures/clip_scatter_re/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/enriched_re/", showWarnings = FALSE, recursive = TRUE)

args = commandArgs(trailingOnly=TRUE)
re_data = read_tsv(args[1])
re_annotations = read_tsv(args[2], col_names = c("chr","start","end","name","score","strand","repeat_name","repeat_class","repeat_family","gc"))  
model_data = read_tsv(args[3])
input_replicate_label = args[4]
clip_replicate_label = args[5]
output_stem = args[6]

logisticb2 = function(x) 1 / (1 + 2**-x)
logitb2 = function(x) log2(x / (1 - x))

re_gc = re_annotations %>% group_by(repeat_name) %>% summarize(gc = mean(gc)) 
model_overdispersion = median(logitb2(model_data$rho))
selected_re_data = re_data %>% rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)) %>% 
	select(-matches("(IP|IN)_[0-9]+$"))	%>% filter(input + clip > 0) %>% inner_join(re_gc) %>% 
	mutate(gc_bin = cut_number(gc,20)) %>% group_by(gc_bin) %>%
	mutate(baseline_l2or = median( logitb2((clip + 1) / (clip + input + 2)) )) %>% ungroup

p_data = selected_re_data %>% group_by(gc_bin, baseline_l2or, clip, input) %>%
	summarize %>% 
	mutate(enrichment_l2or = log2((clip + logisticb2(baseline_l2or)) / (input + 1 - logisticb2(baseline_l2or))) - baseline_l2or) %>%
	mutate(pvalue = pmax(1e-10, 1 - VGAM::pbetabinom(q = clip - 1, size = clip + input, prob = logisticb2(baseline_l2or), rho = model_overdispersion %>% logisticb2 ))) %>%
	inner_join(selected_re_data,.)

q_data = p_data %>% mutate(simple = ifelse(grepl("\\)n",repeat_name), "kmer", "element")) %>% group_by(above_threshold = input + clip >= 10) %>% 
	mutate(qvalue = ifelse(above_threshold, p.adjust(pvalue,"fdr"), NA)) %>% ungroup %>% arrange(pvalue, desc(enrichment_l2or))

label_data = q_data %>% filter(!grepl("\\)n$", repeat_name), enrichment_l2or > 2.5) %>% head(4)
enrichment_max = q_data$enrichment_l2or %>% abs %>% max
pdf(paste0('output/figures/clip_scatter_re/', output_stem, '.clip_test_distribution.pdf'),height = 1.8, width = 2.8)
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

write_tsv(q_data %>% select(-above_threshold), paste0("output/enriched_re/", output_stem, ".enriched_re.tsv.gz"))
