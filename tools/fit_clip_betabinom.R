#!/home/eboyle/bin/Rscript --vanilla  	
# qsub -d $(pwd) -l walltime=35:00 -l nodes=1:ppn=1 -q home-yeo -F 203_HNRNPC_HepG2 fit_clip_betabinom.R
# snakemake -j 20 -s SnakeCLIP_regions.py -k --cluster "qsub -q home-yeo -l walltime={resources.runtime} -N {params.job_name} -e {params.error_out_file} "
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
dir.create("output/figures/clip_distributions/", showWarnings = FALSE, recursive = TRUE)
dir.create("output/clip_model_coef/", showWarnings = FALSE, recursive = TRUE)

count_data = read_tsv(args[1])
experiment = args[2]
given_clip_replicate = args[3]

other_clip_replicates = setdiff(grep("IP_[0-9]*$",names(count_data), value=TRUE), given_clip_replicate)

count_gc_data = count_data[select(count_data, matches("(IP|IN)_[0-9]*$")) %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(gc,10))

clip_betabinom_fit_data = lapply(other_clip_replicates, function(other_clip_replicate)
	{
		clip_gc_data = count_gc_data[as.logical(count_gc_data[,given_clip_replicate] + count_gc_data[,other_clip_replicate] > 0),]
		processed_clip_data = clip_gc_data %>% summarize(given_clip_fraction = mean(.data[[given_clip_replicate]] / (.data[[other_clip_replicate]] + .data[[given_clip_replicate]])) %>% (VGAM::logitlink)) %>% inner_join(clip_gc_data)
		betabinom_fit = VGAM::vglm(cbind(processed_clip_data[[given_clip_replicate]], processed_clip_data[[other_clip_replicate]]) ~ processed_clip_data$given_clip_fraction, VGAM::betabinomial, trace = TRUE)
		betabinom_coefs = betabinom_fit %>% (VGAM::Coef)
		global_given_clip_fraction = sum(clip_gc_data[[given_clip_replicate]]) / sum(clip_gc_data[[other_clip_replicate]] + clip_gc_data[[given_clip_replicate]])
		distribution_data = processed_clip_data %>% group_by(clip_total = .data[[given_clip_replicate]] + .data[[other_clip_replicate]], .data[[given_clip_replicate]], given_clip_fraction, gc_bin) %>%
			count(name = "count") %>% group_by(clip_total, gc_bin) %>% mutate(pdf = count / sum(count)) %>% rowwise %>% 
			mutate(binomial = dbinom(x = .data[[given_clip_replicate]], size = clip_total, prob = given_clip_fraction %>% (VGAM::logitlink)(inverse=TRUE)),
				betabinomial = VGAM::dbetabinom(x = .data[[given_clip_replicate]], size = clip_total, prob = (betabinom_coefs[1] + given_clip_fraction * betabinom_coefs[3]) %>% (VGAM::logitlink)(inverse=TRUE), rho = betabinom_coefs[2] %>% (VGAM::logitlink)(inverse=TRUE) )) 
		pdf(paste0('output/figures/clip_distributions/', experiment, ".", given_clip_replicate, ".", other_clip_replicate, '.clip_distribution.pdf'),height = 3.5, width = 6)
		print(
			ggplot(distribution_data %>% filter(clip_total <=6) %>% 
				pivot_longer(values_to="expected",names_to="distribution",-clip_total:-pdf) %>% 
				mutate(gc_bin = as.numeric(gc_bin) * 10), 
				aes(.data[[given_clip_replicate]] / (clip_total), pdf)) +
			geom_area(data=distribution_data %>% filter(clip_total <=6) %>% mutate(gc_bin = as.numeric(gc_bin) * 10)) +
			theme_bw(base_size = 7)+
			geom_line(aes(y = expected,color = distribution),linetype="21") + 
			facet_grid(clip_total~gc_bin) + xlab(paste0("Fraction of clip reads: ", given_clip_replicate)) + ylab("PDF") + 
			scale_x_continuous(breaks = c(0,1)) + ggtitle("Binned by GC percentile")
		)
		dev.off()

		betabinom_coefs %>% as_tibble(rownames="coef") %>% transmute(coef = c("mu","rho","gc_bias_term"),value) %>% pivot_wider(names_from=coef,values_from=value) %>% mutate(global_given_replicate_fraction = global_given_clip_fraction, given_replicate = given_clip_replicate, other_replicate = other_clip_replicate, experiment = experiment)
	}
) %>% bind_rows

write_tsv(clip_betabinom_fit_data, paste0("output/clip_model_coef/", experiment, ".", given_clip_replicate, ".tsv"))
