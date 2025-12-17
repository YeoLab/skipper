# Load tidyverse for data manipulation and I/O.
library(tidyverse)

# Parse command-line arguments.
args = commandArgs(trailingOnly=TRUE)

# Ensure output directories exist for figures and coefficient tables.
dir.create("output/secondary_figures/figures/clip_distributions/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/secondary_results/clip_model_coef/", showWarnings = FALSE, recursive = TRUE)

# Input data:
# args[1] = count table (with CLIP and input counts, GC bins, etc.)
# args[2] = experiment label
# args[3] = the given CLIP replicate to compare against others
count_data = read_tsv(args[1])
experiment = args[2]
given_clip_replicate = args[3]

# Identify the “other” CLIP replicates (all columns ending in IP_<number> except the given one).
other_clip_replicates = setdiff(grep("IP_[0-9]*$",names(count_data), value=TRUE), given_clip_replicate)

# Subset count data to rows with nonzero total counts across replicates (input+clip > 0),
# then group into GC deciles for modeling GC-dependent biases.
count_gc_data = count_data[select(count_data, matches("(IP|IN)_[0-9]*$")) %>% rowSums > 0,] %>% 
  group_by(gc_bin = cut_number(gc,10))

# Fit beta-binomial models comparing the given replicate vs. each other replicate.
clip_betabinom_fit_data = lapply(other_clip_replicates, function(other_clip_replicate)
{
  # Restrict to rows where this replicate pair has nonzero total counts.
  clip_gc_data = count_gc_data[as.logical(count_gc_data[,given_clip_replicate] + count_gc_data[,other_clip_replicate] > 0),]

  # Calculate the logit of the fraction for the given replicate within the pair, per GC bin.
  processed_clip_data = clip_gc_data %>% 
    summarize(given_clip_fraction = mean(.data[[given_clip_replicate]] / (.data[[other_clip_replicate]] + .data[[given_clip_replicate]])) %>% (VGAM::logitlink)) %>% 
    inner_join(clip_gc_data)

  # Fit a beta-binomial regression with GC-bias term as predictor.
  betabinom_fit = VGAM::vglm(
    cbind(processed_clip_data[[given_clip_replicate]], processed_clip_data[[other_clip_replicate]]) ~ processed_clip_data$given_clip_fraction, 
    VGAM::betabinomial, trace = TRUE
  )

  # Extract coefficients (mu, rho, GC-bias term).
  betabinom_coefs = betabinom_fit %>% (VGAM::Coef)

  # Compute global fraction of reads in the given replicate across this pair.
  global_given_clip_fraction = sum(clip_gc_data[[given_clip_replicate]]) / sum(clip_gc_data[[other_clip_replicate]] + clip_gc_data[[given_clip_replicate]])

  # Build distribution data:
  # Count observed outcomes per total and GC bin, normalize to PDF,
  # then compare against fitted binomial and beta-binomial distributions
  distribution_data = processed_clip_data %>% 
    group_by(clip_total = .data[[given_clip_replicate]] + .data[[other_clip_replicate]], .data[[given_clip_replicate]], given_clip_fraction, gc_bin) %>%
    count(name = "count") %>% 
    group_by(clip_total, gc_bin) %>% mutate(pdf = count / sum(count)) %>% 
    rowwise %>% 
    mutate(
      binomial = dbinom(x = .data[[given_clip_replicate]], size = clip_total, prob = given_clip_fraction %>% (VGAM::logitlink)(inverse=TRUE)),
      betabinomial = VGAM::dbetabinom(
        x = .data[[given_clip_replicate]], 
        size = clip_total, 
        prob = (betabinom_coefs[1] + given_clip_fraction * betabinom_coefs[3]) %>% (VGAM::logitlink)(inverse=TRUE), 
        rho = betabinom_coefs[2] %>% (VGAM::logitlink)(inverse=TRUE)
      )
    )

  # Plot empirical vs. fitted distributions for small totals (≤ 6), binned by GC decile.
  pdf(paste0('output/secondary_figures/figures/clip_distributions/', experiment, ".", given_clip_replicate, ".", other_clip_replicate, '.clip_distribution.pdf'),height = 3.5, width = 6)
  print(
    ggplot(distribution_data %>% filter(clip_total <=6) %>% 
             pivot_longer(values_to="expected",names_to="distribution",-clip_total:-pdf) %>% 
             mutate(gc_bin = as.numeric(gc_bin) * 10), 
           aes(.data[[given_clip_replicate]] / (clip_total), pdf)) +
      geom_area(data=distribution_data %>% filter(clip_total <=6) %>% mutate(gc_bin = as.numeric(gc_bin) * 10)) +
      theme_bw(base_size = 7)+
      geom_line(aes(y = expected,color = distribution),linetype="21") + 
      facet_grid(clip_total~gc_bin) + 
      xlab(paste0("Fraction of clip reads: ", given_clip_replicate)) + ylab("PDF") + 
      scale_x_continuous(breaks = c(0,1)) + ggtitle("Binned by GC percentile")
  )
  dev.off()

  # Save coefficients and bookkeeping for this replicate pair as a tidy tibble.
  betabinom_coefs %>% as_tibble(rownames="coef") %>% 
    transmute(coef = c("mu","rho","gc_bias_term"),value) %>% 
    pivot_wider(names_from=coef,values_from=value) %>% 
    mutate(
      global_given_replicate_fraction = global_given_clip_fraction, 
      given_replicate = given_clip_replicate, 
      other_replicate = other_clip_replicate, 
      experiment = experiment
    )
}
) %>% bind_rows

# Write all fitted coefficients across replicate pairs to file.
write_tsv(clip_betabinom_fit_data, paste0("output/secondary_results/clip_model_coef/", experiment, ".", given_clip_replicate, ".tsv"))
