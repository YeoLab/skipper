library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)

source_mode         <- args[1]
input_gff           <- args[2]
master_ranking_path <- args[3]
output_gff          <- args[4]
output_accession    <- args[5]

# ---------------------- Helpers: I/O & mapping ----------------------

# Read GFF3 into a data.frame. 
read_gff_as_df <- function(path) {
    gr = rtracklayer::readGFF(path)
    gr
}

# Convert back to GRanges and export as GFF3 (attributes = remaining columns)
write_df_as_gff3 <- function(df, path) {
  gr = makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  rtracklayer::export(gr, path, format = "gff3")
}

# process_gencode
process_gencode <- function(df, master_ranking_path, ranking_output, gff_output) {

    # Extract out the transcript support levels.
    tsl <- as.character(df$transcript_support_level)

    # Some transcript support levels have a string that says NA
    is_literal_NA <- !is.na(tsl) & tsl == "NA"

    # Convert tsl to real numbers and convert all NAs to true missing values.  
    tsl_num <- rep(NA_real_, length(tsl))
    suppressWarnings({
        tsl_num[!is_literal_NA] <- as.numeric(tsl[!is_literal_NA])
    })

    # Create a list including everything we want to keep. 
    keep <- (is_literal_NA | is.na(tsl_num) | (tsl_num <= 3))

    # filter the dataframe to only include the contents of keep. 
    df_filtered <- df[keep, , drop = FALSE]
    rownames(df_filtered) <- NULL
    # (Need to keep the NAs as these often represent RNA types outside of the support level structure.) 

    # Find and save all unique transcript types. 
    for_accession <- data.frame(
        accession = unique(df_filtered$transcript_type),
        stringsAsFactors = FALSE
    )

    # Load up the master accession rankings. 
    master <- readr::read_tsv(master_ranking_path, show_col_types = FALSE)

    # Filter the accession rankings to only include those actually in the GFF file. 
    acession_ranking <- left_join(for_accession, master, by = "accession") |>
        tidyr::drop_na() |>
        arrange(rank)

    # Remake the ranking numbers. 
    acession_ranking$rank <- seq_len(nrow(acession_ranking))

    # Save the filtered accession ranking. 
    readr::write_tsv(acession_ranking, ranking_output)

    # Save the filtered gff to file. 
    write_df_as_gff3(df_filtered, gff_output)
}

# Load in the gff file. 
df <- read_gff_as_df(input_gff)

if (tolower(source_mode) == "ensembl") {
  #process_ensembl(df, master_ranking_path, output_accession, output_gff)
} else if (tolower(source_mode) == "gencode") {
  process_gencode(df, master_ranking_path, output_accession, output_gff)
} else {
  stop('Unsupported --source: ', source_mode, ' (expected "ensembl" or "gencode")')
}
