# Load required libraries.
library(tidyverse)
library(ggrepel)

# Print working directory (useful for debugging file paths when running in Snakemake or cluster jobs).
print(getwd())

# Parse command-line arguments.
args = commandArgs(trailingOnly=TRUE)

# Create output directories if they don’t already exist.
dir.create("output/figures/secondary_figures/input_distributions/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/input_model_coef/", showWarnings = FALSE, recursive = TRUE)

# Inputs:
# args[1] = count data table
# args[2] = experiment label
# args[3] = given input replicate name (column)
count_data = read_tsv(args[1])
experiment = args[2]
given_input_replicate = args[3]

# Identify all other input replicates (columns named IN_<number>, excluding the given one).
other_input_replicates = setdiff(grep("IN_[0-9]+$", names(count_data), value=TRUE), given_input_replicate)

# Filter to rows with nonzero reads across all replicates, then group by GC bin (10 equal-frequency bins).
count_gc_data = count_data[select(count_data, matches("(IP|IN)_[0-9]*$")) %>% rowSums > 0,] %>%
  group_by(gc_bin = cut_number(gc,10))

# Fit beta-binomial models for each pair: given input replicate vs. one other input replicate.
input_betabinom_fit_data = lapply(other_input_replicates, function(other_input_replicate) {
    
    # Keep rows where this replicate pair has nonzero total reads.
    input_gc_data = count_gc_data[as.logical(count_gc_data[,given_input_replicate] + count_gc_data[,other_input_replicate] > 0),]
    
    # Compute logit-transformed mean fraction of reads assigned to the given replicate.
    processed_input_data = input_gc_data %>%
      summarize(given_input_fraction = mean(.data[[given_input_replicate]] / 
                                            (.data[[other_input_replicate]] + .data[[given_input_replicate]])) %>% 
                  (VGAM::logitlink)) %>% 
      inner_join(input_gc_data)
    
    # Fit beta-binomial regression model with GC-bias covariate.
    betabinom_fit = VGAM::vglm(
      cbind(processed_input_data[[given_input_replicate]], processed_input_data[[other_input_replicate]]) ~ processed_input_data$given_input_fraction,
      VGAM::betabinomial, trace = TRUE
    )
    betabinom_coefs = betabinom_fit %>% (VGAM::Coef)
    
    # Global fraction across all sites (used for baseline comparisons).
    global_given_input_fraction = sum(input_gc_data[[given_input_replicate]]) / 
                                  sum(input_gc_data[[other_input_replicate]] + input_gc_data[[given_input_replicate]])
    
    # Build distribution of observed counts and compare to expected (binomial vs beta-binomial).
    distribution_data = processed_input_data %>%
      group_by(input_total = .data[[given_input_replicate]] + .data[[other_input_replicate]], 
               .data[[given_input_replicate]], given_input_fraction, gc_bin) %>%
      count(name = "count") %>%
      group_by(input_total, gc_bin) %>%
      mutate(pdf = count / sum(count)) %>% 
      rowwise %>%
      mutate(
        binomial = dbinom(
          x = .data[[given_input_replicate]], 
          size = input_total, 
          prob = given_input_fraction %>% (VGAM::logitlink)(inverse=TRUE)
        ),
        betabinomial = VGAM::dbetabinom(
          x = .data[[given_input_replicate]], 
          size = input_total, 
          prob = (betabinom_coefs[1] + given_input_fraction * betabinom_coefs[3]) %>% (VGAM::logitlink)(inverse=TRUE), 
          rho = betabinom_coefs[2] %>% (VGAM::logitlink)(inverse=TRUE)
        )
      ) 
    
    # Plot observed vs. expected distributions, faceted by GC bin and total read count.
    pdf(paste0("output/figures/secondary_figures/input_distributions/", experiment, ".", given_input_replicate, ".", other_input_replicate, '.input_distribution.pdf'),
        height = 3.5, width = 6)
    print(
      ggplot(distribution_data %>% filter(input_total <= 6) %>% 
               pivot_longer(values_to="expected", names_to="distribution", -input_total:-pdf) %>%
               mutate(gc_bin = as.numeric(gc_bin) * 10), 
             aes(.data[[given_input_replicate]] / (input_total), pdf)) +
        geom_area(data=distribution_data %>% filter(input_total <= 6) %>% mutate(gc_bin = as.numeric(gc_bin) * 10)) +
        theme_bw(base_size = 7) +
        geom_line(aes(y = expected, color = distribution), linetype="21") + 
        facet_grid(input_total ~ gc_bin) + 
        xlab(paste0("Fraction of input reads: ", given_input_replicate)) + 
        ylab("PDF") + 
        scale_x_continuous(breaks = c(0,1)) + 
        ggtitle("Binned by GC percentile")
    )
    dev.off()
    
    # Collect coefficients in tidy format.
    betabinom_coefs %>%
      as_tibble(rownames="coef") %>%
      transmute(coef = c("mu","rho","gc_bias_term"), value) %>%
      pivot_wider(names_from=coef, values_from=value) %>%
      mutate(global_given_replicate_fraction = global_given_input_fraction,
             given_replicate = given_input_replicate,
             other_replicate = other_input_replicate,
             experiment = experiment)
}
) %>% bind_rows

# Save results to file.
write_tsv(input_betabinom_fit_data, paste0("output/secondary_results/input_model_coef/", experiment, ".", given_input_replicate, ".tsv"))
