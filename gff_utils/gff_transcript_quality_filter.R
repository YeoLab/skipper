# Load core libraries for data manipulation/IO
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# Parse command-line arguments in the expected order.
args <- commandArgs(trailingOnly = TRUE)
source_mode         <- args[1]
input_gff           <- args[2]
output_gff          <- args[3]

# Read GFF3 into a data.frame.
read_gff_as_df <- function(path) {
    df <- rtracklayer::readGFF(path)

    # Make sure we are working with a plain data.frame.
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    df
}

# Convert back to GRanges and export as GFF3.
write_df_as_gff3 <- function(df, path) {

    # Drop row names.
    rownames(df) <- NULL

    # Remove helper columns if they somehow made it through.
    helper_cols <- c(
        "tag_str",
        "is_mane",
        "is_canonical",
        "is_primary",
        "is_appris",
        "is_ccds",
        "is_basic",
        "priority",
        "tie_break_score",
        "min_priority",
        "max_tie_break_score"
    )
    df <- df[, setdiff(colnames(df), helper_cols), drop = FALSE]

    # Coerce list columns to semicolon-delimited character strings.
    df[] <- lapply(df, function(col) {
        if (is.list(col)) {
            vapply(col, function(x) {
                if (length(x) == 0 || all(is.na(x))) {
                    NA_character_
                } else {
                    paste(as.character(x), collapse = ";")
                }
            }, character(1))
        } else {
            col
        }
    })

    # Ensure core GFF columns exist with expected names/types.
    required_cols <- c("seqid", "start", "end")
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop(
            paste0(
                "Missing required columns for GFF export: ",
                paste(missing_cols, collapse = ", ")
            )
        )
    }

    if (!"strand" %in% colnames(df)) {
        df$strand <- "*"
    }
    if (!"phase" %in% colnames(df)) {
        df$phase <- NA
    }

    df$start <- as.integer(df$start)
    df$end <- as.integer(df$end)
    df$strand <- as.character(df$strand)
    df$phase <- suppressWarnings(as.integer(df$phase))

    gr <- makeGRangesFromDataFrame(
        df,
        seqnames.field = "seqid",
        start.field = "start",
        end.field = "end",
        strand.field = "strand",
        keep.extra.columns = TRUE,
        ignore.strand = FALSE
    )

    con <- gzfile(path, "w")
    rtracklayer::export(gr, con, format = "gff3")
    close(con)
}

# process_gencode
process_gencode <- function(df, gff_output) {

    # Extract transcript support levels.
    tsl <- as.character(df$transcript_support_level)

    # Some transcript support levels are literally the string "NA".
    is_literal_NA <- !is.na(tsl) & tsl == "NA"

    # Convert tsl to numeric where possible.
    tsl_num <- rep(NA_real_, length(tsl))
    suppressWarnings({
        tsl_num[!is_literal_NA] <- as.numeric(tsl[!is_literal_NA])
    })

    # Keep transcripts with TSL <= 3, plus true/literal NAs.
    keep <- is_literal_NA | is.na(tsl_num) | (tsl_num <= 3)

    df_filtered <- df[keep, , drop = FALSE]
    rownames(df_filtered) <- NULL

    # Build transcript-level priority table.
    df_tx <- df_filtered %>%
        filter(type == "transcript") %>%
        mutate(
            tag_str = vapply(tag, function(x) {
                if (length(x) == 0 || all(is.na(x))) {
                    ""
                } else {
                    paste(as.character(x), collapse = ";")
                }
            }, character(1)),

            is_mane = str_detect(tag_str, "MANE_Select"),
            is_canonical = str_detect(tag_str, "Ensembl_canonical"),
            is_primary = str_detect(tag_str, "GENCODE_Primary"),
            is_appris = str_detect(tag_str, "appris_principal_[1-5]"),
            is_ccds = str_detect(tag_str, "CCDS"),
            is_basic = str_detect(tag_str, "(^|;)basic($|;)"),

            priority = case_when(
                is_mane ~ 1L,
                is_canonical ~ 2L,
                is_primary ~ 3L,
                is_appris ~ 4L,
                is_ccds ~ 5L,
                is_basic ~ 6L,
                TRUE ~ 7L
            ),

            tie_break_score = case_when(
                priority == 1L ~ 0L,
                priority == 2L ~ as.integer(is_primary) + as.integer(is_appris) + as.integer(is_ccds) + as.integer(is_basic),
                priority == 3L ~ as.integer(is_appris) + as.integer(is_ccds) + as.integer(is_basic),
                priority == 4L ~ as.integer(is_ccds) + as.integer(is_basic),
                priority == 5L ~ as.integer(is_basic),
                TRUE ~ 0L
            )
        ) %>%
        group_by(gene_id) %>%
        mutate(min_priority = min(priority, na.rm = TRUE)) %>%
        ungroup()

    # Keep all top transcripts after tie-breaking.
    best_tx <- df_tx %>%
        filter(priority == min_priority) %>%
        group_by(gene_id) %>%
        mutate(max_tie_break_score = max(tie_break_score, na.rm = TRUE)) %>%
        filter(tie_break_score == max_tie_break_score) %>%
        ungroup()

    # Get selected transcript IDs.
    selected_tx_ids <- unique(best_tx$transcript_id)

    # Keep all rows belonging to selected transcripts, plus gene rows for selected genes.
    selected_gene_ids <- unique(best_tx$gene_id)

    df_out <- df_filtered %>%
        filter(
            gene_id %in% selected_gene_ids &
            (
                type == "gene" |
                transcript_id %in% selected_tx_ids |
                (type == "transcript" & transcript_id %in% selected_tx_ids)
            )
        )

    # Write filtered GFF.
    write_df_as_gff3(df_out, gff_output)
}

# Load in the gff file. 
df <- read_gff_as_df(input_gff)

if (tolower(source_mode) == "ensembl") {
  df <- df |> as.data.frame()
  process_ensembl(df, output_gff)
} else if (tolower(source_mode) == "gencode") {
  process_gencode(df, output_gff)
} else {
  stop('Unsupported --source: ', source_mode, ' (expected "ensembl" or "gencode")')
}