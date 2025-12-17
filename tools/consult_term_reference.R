# Load core libraries for data manipulation/IO (tidyverse), plot layout (cowplot), dendrogram handling (ggdendro), and color scales (viridis).
library(tidyverse)
library(cowplot)
library(ggdendro)
library(viridis)

# Parse command-line arguments in the expected order.
args = commandArgs(trailingOnly=TRUE)
enriched_window_file = args[1] 
term_gmt_file = args[2] 
term_reference_file = args[3]
jaccard_rds_file = args[4]
prefix = args[5]

# Ensure output directories exist for tables and figures.
dir.create("output/gene_sets/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/gene_sets/", showWarnings = FALSE, recursive = TRUE)

# Load input tables and gene sets.
enriched_windows = read_tsv(enriched_window_file)
terms = fgsea::gmtPathways(term_gmt_file)
term_reference = read_tsv(term_reference_file)

# Keep only those terms that appear in the reference table to align with available priors.
valid_terms = terms[unique(term_reference$term)]

# Count enriched windows per gene (splitting multi-gene entries on “:”).
# Produces a tibble with columns: gene, n_windows_enriched.
clip_gene_data = enriched_windows %>% summarize(table(gene = strsplit(gene_name,split=":") %>% unlist) %>% as_tibble) %>% rename(n_windows_enriched = n)

# Total number of enriched windows across all genes, used as the binomial denominator.
n_windows_total = clip_gene_data %>% with(sum(n_windows_enriched))

# Compute the size of each valid term (number of genes in the set).
term_sizes = lapply(valid_terms, function(term_list) length(term_list)) %>% as_tibble %>% pivot_longer(all_of(names(.)),names_to="term",values_to="n_genes_term") 

# For each term, count how many genes are enriched and how many enriched windows those genes account for.
# Join sizes and replace missing counts with zeros for terms that have zero overlap with enriched genes.
term_data = lapply(valid_terms, function(term_list) filter(clip_gene_data, gene %in% term_list) %>% summarize(n_genes_enriched = n(), n_windows_enriched = sum(n_windows_enriched))) %>% bind_rows(.id = "term") %>% 
	left_join(term_sizes,.) %>% 
	mutate(n_genes_enriched = ifelse(is.na(n_genes_enriched), 0, n_genes_enriched),n_windows_enriched = ifelse(is.na(n_windows_enriched), 0, n_windows_enriched)) 

# Perform a binomial test per term comparing observed windows to expectation (from reference fraction f_term_windows).
# Compute a simple FDR (Bonferroni-like) by p * number_of_terms, clamped to 1, and add log2 fold-change of observed vs expected.
binom_data = term_data %>% inner_join(term_reference) %>% rowwise %>% 
	summarize((broom::tidy)(binom.test(x= n_windows_enriched, n = n_windows_total, p = f_term_windows)) %>% ungroup %>% mutate(term=term, l2fc = log2(estimate/f_term_windows))) %>% 
	inner_join(term_data %>% mutate(n_windows_total = n_windows_total),.) %>% mutate(p_unadjusted = ifelse(p.value == 0, min(p.value[p.value!=0]), p.value), p_adj = pmin(1,p_unadjusted * nrow(.))) %>% select(-p.value) %>% arrange(p_unadjusted,l2fc %>% desc)

# Save the term-level enrichment results.
write_tsv(binom_data, paste0("output/gene_sets/", prefix, ".enriched_terms.tsv.gz"))

# Load precomputed term–term Jaccard similarity matrix for clustering and visualization of the top terms.
jaccard_mat = readRDS(jaccard_rds_file)

# Choose the top significant terms for plotting (positive L2FC and FDR < 0.05, capped at 15).
top_terms = binom_data %>% filter(l2fc > 0,p_adj < 0.05) %>% head(15) %>% pull(term)

# If we have fewer than two terms, emit a placeholder figure and exit early.
if(length(top_terms) < 2) {
	pdf(paste0("output/figures/gene_sets/", prefix, ".clustered_top_terms.pdf"),height= 1,width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "0 or 1 significant terms") + theme_void())
	dev.off()
	quit()
}

# Hierarchical clustering of the selected terms using the Jaccard distance derived from the similarity matrix.
top_terms_hc = jaccard_mat[top_terms,top_terms] %>% as.dist %>% hclust
top_terms_hcd = as.dendrogram(top_terms_hc)
top_terms_hcd_data = dendro_data(top_terms_hcd, type = "rectangle")

# Bar plot showing L2FC per term, ordered by the dendrogram’s label order, with -log10(FDR) as fill and gene counts as labels.
bars = binom_data %>% filter(term %in% top_terms) %>% 
	mutate(term = factor(term, levels = top_terms_hc$labels)) %>%
	ggplot(aes(x=term,y=l2fc, fill = -log10(p_adj), label = paste0(n_genes_enriched, "/", n_genes_term))) + theme_bw(base_size = 6) + 
	theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.margin=margin(5,5,5,0)) +
	geom_bar(stat="identity", width=0.9) + geom_text(size = 2, color = "white", hjust = 1.1) + coord_flip() + scale_x_discrete(position = "top") + 
	xlab("") + ylab("L2FC O/E windows") + 
	scale_fill_viridis(option="B",end = 0.92)

# Dendrogram panel for the same set of terms to display their similarity structure.
phylogeny = ggplot(segment(top_terms_hcd_data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) + theme_dendro() + theme(plot.margin=margin(5,0,5,5))

# Combine dendrogram and bars side by side and write to PDF.
pdf(paste0("output/figures/gene_sets/", prefix, ".clustered_top_terms.pdf"),height= 1.9,width = 6)
plot_grid(phylogeny,bars,rel_widths=c(.25,1))
dev.off()

