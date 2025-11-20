# Load core tidyverse for data manipulation and I/O, and ggrepel for non-overlapping labels in plots.
library(tidyverse)
library(ggrepel)

# Parse CLI args and ensure output directories exist for figures and coefficient tables.
args = commandArgs(trailingOnly=TRUE)
dir.create("output/figures/secondary_figures/clip_scatter_re/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/clip_model_coef_re/", showWarnings = FALSE, recursive = TRUE)

# Inputs: repeat-element count table, experiment label, and the “given” CLIP replicate to model against others.
re_data = read_tsv(args[1])
experiment = args[2]
given_clip_replicate = args[3]

# Identify the “other” CLIP replicates to compare against (columns named like IP_<number> except the given one).
other_clip_replicates = setdiff(grep("IP_[0-9]*$",names(re_data), value=TRUE), given_clip_replicate)

# Fit a beta-binomial null per pair (given vs. each other replicate), generate a diagnostic plot, and collect coefficients.
clip_betabinom_re_fit_data = lapply(other_clip_replicates, function(other_clip_replicate)
    {
        # Restrict to rows with nonzero total CLIP counts across the pair to avoid degenerate sizes.
        nonzero_re_data = re_data[as.logical(re_data[,given_clip_replicate] + re_data[,other_clip_replicate] > 0),]

        # Empirical fraction of reads in the given replicate across the pair, used for visualization of enrichment.
        given_clip_fraction = sum(nonzero_re_data[[given_clip_replicate]]) / sum(nonzero_re_data[[other_clip_replicate]] + nonzero_re_data[[given_clip_replicate]])

        # Fit beta-binomial with intercept only (global mu and rho) using VGAM’s vglm + betabinomial family.
        betabinom_fit = VGAM::vglm(cbind(nonzero_re_data[[given_clip_replicate]], nonzero_re_data[[other_clip_replicate]]) ~ 1, VGAM::betabinomial, trace = TRUE)

        # Extract coefficients: “mu” (mean success probability) and “rho” (intra-class correlation / overdispersion).
        betabinom_coefs = betabinom_fit %>% (VGAM::Coef)

        # Collapse to unique count pairs and compute null-model enrichment and p-values under the fitted beta-binomial.
        p_data_null = nonzero_re_data %>% group_by(.data[[given_clip_replicate]], .data[[other_clip_replicate]]) %>%
            summarize %>% mutate(enrichment_l2or = log2((.data[[given_clip_replicate]] + betabinom_coefs["mu"]) / (.data[[other_clip_replicate]] + 1 - betabinom_coefs["mu"])) - log2(betabinom_coefs["mu"] / (1 - betabinom_coefs["mu"]))) %>%
            mutate(pvalue = pmax(1e-13, 1 - VGAM::pbetabinom(q = .data[[given_clip_replicate]] - 1, size = .data[[given_clip_replicate]] + .data[[other_clip_replicate]], prob = betabinom_coefs["mu"], rho = betabinom_coefs["rho"]))) %>%
            inner_join(nonzero_re_data,.) 

        # FDR within the subset with adequate total counts (≥ 10) to stabilize tail probabilities.
        q_data_null = p_data_null[as.logical(p_data_null[,given_clip_replicate] + p_data_null[,other_clip_replicate] >= 10),] %>% mutate(qvalue = p.adjust(pvalue, "fdr"))

        # Pick a few top non-kmer labels for annotation in the diagnostic plot.
        label_data = q_data_null %>% arrange(qvalue) %>% filter(!grepl("\\)n$", repeat_name)) %>% head(4)

        # Diagnostic scatter: total counts vs. log2 enrichment, colored by q-value threshold, with the null l2(μ/(1−μ)) reference line.
        pdf(paste0('output/figures/secondary_figures/clip_scatter_re/', experiment, ".", given_clip_replicate, ".", other_clip_replicate, '.clip_null_distribution.pdf'),height = 1.8, width = 2.8)
        print(
            nonzero_re_data %>% inner_join(q_data_null) %>% 
            ggplot(aes(.data[[given_clip_replicate]]+.data[[other_clip_replicate]], log2((given_clip_fraction+.data[[given_clip_replicate]]) / (1 - given_clip_fraction +.data[[other_clip_replicate]])) - log2(given_clip_fraction / (1 - given_clip_fraction)), color = ifelse(qvalue < 0.05, "q < 0.05", "NS"))) + 
                geom_hline(yintercept = log2(betabinom_coefs["mu"] / (1 - betabinom_coefs["mu"])) ) +
                theme_bw(base_size = 7) +
                theme(legend.title=element_blank()) +
                geom_point(size = 0.8, stroke = 0) + scale_x_log10() + xlab("Total reads") + ylab("Log2 enrichment") + 
                geom_text_repel(data = label_data, aes(label = repeat_name),min.segment.length = 0.1, color = "black",size = 2) + 
                theme(panel.grid.minor = element_blank()) + scale_color_manual(values = c("#bdbdbd","#fec332"))
            
        )
        dev.off()

        # Return coefficients and bookkeeping per comparison as a tidy row for downstream use.
        betabinom_coefs %>% as_tibble(rownames="coef") %>% transmute(coef = c("mu","rho"),value) %>% pivot_wider(names_from=coef,values_from=value) %>% mutate(given_clip_fraction = given_clip_fraction, given_replicate = given_clip_replicate, other_replicate = other_clip_replicate, experiment = experiment)
    }
) %>% bind_rows

# Write out the per-pair fitted parameters (mu, rho) and the observed fraction to a TSV for reuse in testing steps.
write_tsv(clip_betabinom_re_fit_data, paste0("output/clip_model_coef_re/", experiment, ".", given_clip_replicate, ".tsv"))
