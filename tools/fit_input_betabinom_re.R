library(tidyverse)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)
dir.create("output/figures/input_scatter_re/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/input_model_coef_re/", showWarnings = FALSE, recursive = TRUE)

re_data = read_tsv(args[1])
experiment = args[2]
given_input_replicate = args[3]

other_input_replicates = setdiff(grep("IN_[0-9]*$",names(re_data), value=TRUE), given_input_replicate)

input_betabinom_re_fit_data = lapply(other_input_replicates, function(other_input_replicate)
    {
        nonzero_re_data = re_data[as.logical(re_data[,given_input_replicate] + re_data[,other_input_replicate] > 0),]
        given_input_fraction = sum(nonzero_re_data[[given_input_replicate]]) / sum(nonzero_re_data[[other_input_replicate]] + nonzero_re_data[[given_input_replicate]])
        betabinom_fit = VGAM::vglm(cbind(nonzero_re_data[[given_input_replicate]], nonzero_re_data[[other_input_replicate]]) ~ 1, VGAM::betabinomial, trace = TRUE)
        betabinom_coefs = betabinom_fit %>% (VGAM::Coef)
        p_data_null = nonzero_re_data %>% group_by(.data[[given_input_replicate]], .data[[other_input_replicate]]) %>%
            summarize %>% mutate(enrichment_l2or = log2((.data[[given_input_replicate]] + betabinom_coefs["mu"]) / (.data[[other_input_replicate]] + 1 - betabinom_coefs["mu"])) - log2(betabinom_coefs["mu"] / (1 - betabinom_coefs["mu"]))) %>%
            mutate(pvalue = pmax(1e-13, 1 - VGAM::pbetabinom(q = .data[[given_input_replicate]] - 1, size = .data[[given_input_replicate]] + .data[[other_input_replicate]], prob = betabinom_coefs["mu"], rho = betabinom_coefs["rho"]))) %>%
            inner_join(nonzero_re_data,.) 

        q_data_null = p_data_null[as.logical(p_data_null[,given_input_replicate] + p_data_null[,other_input_replicate] >= 10),] %>% mutate(qvalue = p.adjust(pvalue, "fdr"))

        label_data = q_data_null %>% arrange(qvalue) %>% filter(!grepl("\\)n$", repeat_name)) %>% head(4)
        pdf(paste0('output/figures/input_scatter_re/', experiment, ".", given_input_replicate, ".", other_input_replicate, '.input_null_distribution.pdf'),height = 1.8, width = 2.8)
        print(
            nonzero_re_data %>% inner_join(q_data_null) %>% 
            ggplot(aes(.data[[given_input_replicate]]+.data[[other_input_replicate]], log2((given_input_fraction+.data[[given_input_replicate]]) / (1 - given_input_fraction +.data[[other_input_replicate]])) - log2(given_input_fraction / (1 - given_input_fraction)), color = ifelse(qvalue < 0.05, "q < 0.05", "NS"))) + 
                geom_hline(yintercept = log2(betabinom_coefs["mu"] / (1 - betabinom_coefs["mu"])) ) +
                theme_bw(base_size = 7) +
                theme(legend.title=element_blank()) +
                geom_point(size = 0.8, stroke = 0) + scale_x_log10() + xlab("Total reads") + ylab("Log2 enrichment") + 
                geom_text_repel(data = label_data, aes(label = repeat_name),min.segment.length = 0.1, color = "black",size = 2) + 
                theme(panel.grid.minor = element_blank()) + scale_color_manual(values = c("#bdbdbd","#fec332"))
            
        )
        dev.off()

        betabinom_coefs %>% as_tibble(rownames="coef") %>% transmute(coef = c("mu","rho"),value) %>% pivot_wider(names_from=coef,values_from=value) %>% mutate(given_input_fraction = given_input_fraction, given_replicate = given_input_replicate, other_replicate = other_input_replicate, experiment = experiment)
    }
) %>% bind_rows

write_tsv(input_betabinom_re_fit_data, paste0("output/input_model_coef_re/", experiment, ".", given_input_replicate, ".tsv"))
