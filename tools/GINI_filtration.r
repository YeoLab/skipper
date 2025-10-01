# Load necessary packages. 
library(glue)
library(tidyverse)

# Define the path and condition names for the data. 
args = commandArgs(trailingOnly=TRUE)
input_window_path = args[1]
input_nt_path = args[2]
output_window_path = args[3]

# Load in the window and nucleotide counts data. 
window_data = read_tsv(input_window_path)
nt_counts = read_tsv(nt_path, col_names = c("chr","start","end","name","score","strand","window_n","input","clip"), col_types = c("ciiciciii"))

# Calculate the gini-coefficient of each window from the nt counts. 
window_metrics <- nt_counts %>% 
    group_by(window_n) %>% 
    summarize(
    gini     = DescTools::Gini(clip),
    .groups  = "drop"
    )
    
# add the ratios to nt_counts. 
nt_counts = left_join(nt_counts, window_metrics, by = "window_n")

# Create ID columns for both nt_counts and window_data for merging. 
window_data$ID = paste(window_data$chr, window_data$start, sep = "_")
nt_counts$ID = paste(nt_counts$chr, nt_counts$start, sep = "_")

# Merge the ratio information into the window data
window_data = left_join(window_data, dplyr::select(nt_counts, c("ID", "window_n", "gini")), by = "ID")

# Use the ratios to filter the data. 
window_data_filtered <- window_data %>%
  filter(gini < 0.9)

# Save the 
write_tsv(window_data_filtered, output_window_path)