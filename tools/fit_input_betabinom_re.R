# Load required libraries for data manipulation and plotting.
library(tidyverse)
library(ggrepel)

# Parse command-line arguments.
args = commandArgs(trailingOnly=TRUE)

# Create output directories for figures and model coefficients (if not already present).
dir.create("output/figures/secondary_figures/input_scatter_re/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/input_model_coef_re/", showWarnings = FALSE, recursive = TRUE)

# Inputs:
# args[1] = repeat element (RE) count data
# args[2] = experiment label
# args[3] = the given input replicate to compare against others
re_data = read_tsv(args[1])
experiment = args[2]
given_input_replicate = args[3]

# Identify all other input replicates (columns named IN_<number>, excluding the given one).
other_input_replicates = setdiff(grep("IN_[0-9]*$",names(re_data), value=TRUE), given_input_replicate)

# For each other replicate, fit a beta-binomial model comparing counts against the given input replicate.
input_betabinom_re_fit_data = lapply(other_input_replicates, function(other_input_replicate)
{
    # Keep only rows where the pair of replicates has >0 total reads.
    nonzero_re_data = re_data[as.logical(re_data[,given_input_replicate] + re_data[,other_input_replicate] > 0),]

    # Global fraction of reads in the given replicate across this replicate pair.
    given_input_fraction = sum(nonzero_re_data[[given_input_replicate]]) / 
                           sum(nonzero_re_data[[other_input_replicate]] + nonzero_re_data[[given_input_replicate]])

    # Fit beta-binomial model using VGAM.
    betabinom_fit = VGAM::vglm(
        cbind(nonzero_re_data[[given_input_replicate]], nonzero_re_data[[other_input_replicate]]) ~ 1,
        VGAM::betabinomial, trace = TRUE
    )
    betabinom_coefs = betabinom_fit %>% (VGAM::Coef)

    # Calculate enrichment (log2 odds ratio) and p-values under null model.
    p_data_null = nonzero_re_data %>% 
        group_by(.data[[given_input_replicate]], .data[[other_input_replicate]]) %>%
        summarize %>% 
        mutate(enrichment_l2or = log2((.data[[given_input_replicate]] + betabinom_coefs["mu"]) / 
                                      (.data[[other_input_replicate]] + 1 - betabinom_coefs["mu"])) - 
                                      log2(betabinom_coefs["mu"] / (1 - betabinom_coefs["mu"]))) %>%
        mutate(pvalue = pmax(
            1e-13, 
            1 - VGAM::pbetabinom(
                q = .data[[given_input_replicate]] - 1, 
                size = .data[[given_input_replicate]] + .data[[other_input_replicate]], 
                prob = betabinom_coefs["mu"], 
                rho = betabinom_coefs["rho"]
            )
        )) %>%
        inner_join(nonzero_re_data,.) 

    # Adjust p-values (FDR) for windows with >= 10 total reads.
    q_data_null = p_data_null[as.logical(p_data_null[,given_input_replicate] + p_data_null[,other_input_replicate] >= 10),] %>% 
        mutate(qvalue = p.adjust(pvalue, "fdr"))

    # Pick top 4 enriched repeat elements to label in plot.
    label_data = q_data_null %>% arrange(qvalue) %>% filter(!grepl("\\)n$", repeat_name)) %>% head(4)

    # Generate scatter plot of enrichment vs total reads.
    pdf(paste0('output/figures/secondary_figures/input_scatter_re/', experiment, ".", given_input_replicate, ".", other_input_replicate, '.input_null_distribution.pdf'),
        height = 1.8, width = 2.8)
    print(
        nonzero_re_data %>% inner_join(q_data_null) %>% 
        ggplot(aes(.data[[given_input_replicate]] + .data[[other_input_replicate]], 
                   log2((given_input_fraction + .data[[given_input_replicate]]) / 
                        (1 - given_input_fraction + .data[[other_input_replicate]])) - 
                   log2(given_input_fraction / (1 - given_input_fraction)), 
                   color = ifelse(qvalue < 0.05, "q < 0.05", "NS"))) + 
            geom_hline(yintercept = log2(betabinom_coefs["mu"] / (1 - betabinom_coefs["mu"])) ) +
            theme_bw(base_size = 7) +
            theme(legend.title=element_blank()) +
            geom_point(size = 0.8, stroke = 0) + scale_x_log10() + 
            xlab("Total reads") + ylab("Log2 enrichment") + 
            geom_text_repel(data = label_data, aes(label = repeat_name),
                            min.segment.length = 0.1, color = "black", size = 2) + 
            theme(panel.grid.minor = element_blank()) + 
            scale_color_manual(values = c("#bdbdbd","#fec332"))
    )
    dev.off()

    # Return coefficients and metadata as a tidy tibble.
    betabinom_coefs %>% as_tibble(rownames="coef") %>% 
        transmute(coef = c("mu","rho"),value) %>% 
        pivot_wider(names_from=coef,values_from=value) %>% 
        mutate(given_input_fraction = given_input_fraction, 
               given_replicate = given_input_replicate, 
               other_replicate = other_input_replicate, 
               experiment = experiment)
}
) %>% bind_rows

# Save all fitted coefficient results to file.
write_tsv(input_betabinom_re_fit_data, paste0("output/secondary_results/input_model_coef_re/", experiment, ".", given_input_replicate, ".tsv"))
