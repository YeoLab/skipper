library(tidyverse)
library(cowplot)
library(ggdendro)
library(viridis)

args = commandArgs(trailingOnly=TRUE)

enriched_window_file = args[1] 
term_gmt_file = args[2] 
term_reference_file = args[3]
jaccard_rds_file = args[4]
prefix = args[5] # experiment_label

dir.create("output/gene_sets/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/gene_sets/", showWarnings = FALSE, recursive = TRUE)

enriched_windows = read_tsv(enriched_window_file)
terms = fgsea::gmtPathways(term_gmt_file)
term_reference = read_tsv(term_reference_file)

valid_terms = terms[unique(term_reference$term)]

clip_gene_data = enriched_windows %>% summarize(table(gene = strsplit(gene_name,split=":") %>% unlist) %>% as_tibble) %>% rename(n_windows_enriched = n)
n_windows_total = clip_gene_data %>% with(sum(n_windows_enriched))

term_sizes = lapply(valid_terms, function(term_list) length(term_list)) %>% as_tibble %>% pivot_longer(all_of(names(.)),names_to="term",values_to="n_genes_term") 
term_data = lapply(valid_terms, function(term_list) filter(clip_gene_data, gene %in% term_list) %>% summarize(n_genes_enriched = n(), n_windows_enriched = sum(n_windows_enriched))) %>% bind_rows(.id = "term") %>% 
	left_join(term_sizes,.) %>% 
	mutate(n_genes_enriched = ifelse(is.na(n_genes_enriched), 0, n_genes_enriched),n_windows_enriched = ifelse(is.na(n_windows_enriched), 0, n_windows_enriched)) 

binom_data = term_data %>% inner_join(term_reference) %>% rowwise %>% 
	summarize((broom::tidy)(binom.test(x= n_windows_enriched, n = n_windows_total, p = f_term_windows)) %>% ungroup %>% mutate(term=term, l2fc = log2(estimate/f_term_windows))) %>% 
	inner_join(term_data %>% mutate(n_windows_total = n_windows_total),.) %>% mutate(p_unadjusted = ifelse(p.value == 0, min(p.value[p.value!=0]), p.value), p_adj = pmin(1,p_unadjusted * nrow(.))) %>% select(-p.value) %>% arrange(p_unadjusted,l2fc %>% desc)

write_tsv(binom_data, paste0("output/gene_sets/", prefix, ".enriched_terms.tsv.gz"))

jaccard_mat = readRDS(jaccard_rds_file)

top_terms = binom_data %>% filter(l2fc > 0,p_adj < 0.05) %>% head(15) %>% pull(term)
if(length(top_terms) < 2) {
	pdf(paste0("output/figures/gene_sets/", prefix, ".clustered_top_terms.pdf"),height= 1,width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "0 or 1 significant terms") + theme_void())
	dev.off()
	quit()
}
top_terms_hc = jaccard_mat[top_terms,top_terms] %>% as.dist %>% hclust
top_terms_hcd = as.dendrogram(top_terms_hc)
top_terms_hcd_data = dendro_data(top_terms_hcd, type = "rectangle")

bars = binom_data %>% filter(term %in% top_terms) %>% 
	mutate(term = factor(term, levels = top_terms_hc$labels)) %>%
	ggplot(aes(x=term,y=l2fc, fill = -log10(p_adj), label = paste0(n_genes_enriched, "/", n_genes_term))) + theme_bw(base_size = 6) + 
	theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.margin=margin(5,5,5,0)) +
	geom_bar(stat="identity", width=0.9) + geom_text(size = 2, color = "white", hjust = 1.1) + coord_flip() + scale_x_discrete(position = "top") + 
	xlab("") + ylab("L2FC O/E windows") + 
	scale_fill_viridis(option="B",end = 0.92)

phylogeny = ggplot(segment(top_terms_hcd_data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) + theme_dendro() + theme(plot.margin=margin(5,0,5,5))

pdf(paste0("output/figures/gene_sets/", prefix, ".clustered_top_terms.pdf"),height= 1.9,width = 6)
plot_grid(phylogeny,bars,rel_widths=c(.25,1))
dev.off()

