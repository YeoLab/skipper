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

# Create helper function for adding in stop and start codon information (needed for ensembl). 
add_stop_start <- function(cds, to_keep, strand = "+", start = TRUE) {
  cds_strand <- cds %>% filter(strand == strand)
  if (nrow(cds_strand) == 0) return(cds_strand[0, ])

  if (start && strand == "+") {
    df <- cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(start = min(start), .groups = "drop") %>%
      mutate(end = start + 2L, type = "start_codon")
  } else if (!start && strand == "+") {
    df <- cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(end = max(end), .groups = "drop") %>%
      mutate(start = end - 2L, type = "stop_codon")
  } else if (!start && strand == "-") {
    df <- cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(start = min(start), .groups = "drop") %>%
      mutate(end = start + 2L, type = "stop_codon")
  } else if (start && strand == "-") {
    df <- cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(end = max(end), .groups = "drop") %>%
      mutate(start = end - 2L, type = "start_codon")
  } else {
    stop("Unhandled case in add_stop_start")
  }

  df <- df %>%
    mutate(phase = 0L) %>%
    # Reorder columns to match cds_strand (Python kept same columns)
    select(any_of(colnames(cds_strand)))
  df
}

process_ensembl <- function(df, master_ranking_path, ranking_output, gff_output) {
    colnames(df) <- tolower(colnames(df))

    # Define the types.
    type_1 <- c("gene", "ncRNA_gene", "pseudogene")
    type_2 <- c("mRNA", "lnc_RNA", "ncRNA", "miRNA", "snRNA", "snoRNA", "rRNA",
                "pseudogenic_transcript", "V_gene_segment", "scRNA", "C_gene_segment",
                "Y_RNA", "J_gene_segment", "D_gene_segment", "transcript", "unconfirmed_transcript")
    type_3 <- c("exon", "CDS", "five_prime_UTR", "three_prime_UTR")

    # Build mappings.
    gid_to_gtype <- df %>%
        select(gene_id, biotype) %>%
            filter(!is.na(gene_id), !is.na(biotype)) %>%
                distinct()
    gid_to_gtype_vec <- setNames(gid_to_gtype$biotype, gid_to_gtype$gene_id)

    tid_to_gid_ttype <- df %>%
        select(transcript_id, parent, biotype) %>%
            filter(!is.na(transcript_id), !is.na(parent), !is.na(biotype)) %>%
                distinct()
    
    # Subset to include only the Ensembl annotations (important stuff).
    ensg <- df %>% filter(source == "ensembl") %>% select(
        seqid, source, type, start, end, score, strand, phase,
        id, biotype, gene_id, parent, transcript_id, transcript_support_level, 
    )

    # Step 1: Clean transcript_support_level
    ensg$transcript_support_level <- sub("\\s*\\(.*\\)$", "", as.character(ensg$transcript_support_level))
    
    tsl <- as.character(ensg$transcript_support_level)
    is_literal_NA <- !is.na(tsl) & tsl == "NA"
    
    tsl_num <- rep(NA_real_, length(tsl))
    suppressWarnings({
        tsl_num[!is_literal_NA] <- as.numeric(tsl[!is_literal_NA])
    })
    
    # Step 2: Identify bad transcript IDs (no prefix in ID column)
    is_transcript <- ensg$Feature == "transcript"
    keep_transcript <- (is_literal_NA | is.na(tsl_num) | (tsl_num <= 3))
    bad_transcripts <- ensg$transcript_id[!keep_transcript]
    
    # Step 3: Build the parent-style strings
    bad_transcript_parents <- paste0("transcript:", bad_transcripts)
    
    # Step 4: Drop bad transcripts AND their children (use 'parent' column)
    ensg <- ensg[!(ensg$transcript_id %in% bad_transcripts |
                            ensg$parent %in% bad_transcript_parents), , drop = FALSE]
    rownames(ensg) <- NULL
    
    # Pre-allocate columns.
    ensg <- ensg %>%
        mutate(
          gene_type = NA_character_,
          transcript_type = NA_character_,
          gene_id = as.character(gene_id),
          transcript_id = as.character(transcript_id)
        )

    # Gene-level rows.
    mask1 <- ensg$type %in% type_1
    ensg$gene_type[mask1] <- ensg$biotype[mask1]

    # Transcript rows with parent=gene:*
    mask2 <- ensg$type %in% type_2
    if (any(mask2)) {
        gene_ids <- str_remove(ensg$parent[mask2], "^gene:")
        ensg$gene_id[mask2] <- gene_ids
        ensg$gene_type[mask2] <- unname(gid_to_gtype_vec[gene_ids])
        ensg$transcript_type[mask2] <- ensg$biotype[mask2]
    }

    # type rows with parent=transcript:*
    mask3 <- ensg$type %in% type_3
    if (any(mask3)) {
        transcript_ids <- str_remove(ensg$parent[mask3], "^transcript:")
        ensg$transcript_id[mask3] <- transcript_ids
        
        # map transcript_id -> (parent, biotype)
        idx <- match(transcript_ids, tid_to_gid_ttype$transcript_id)
        parent_vec <- tid_to_gid_ttype$parent[idx]
        ttype_vec  <- tid_to_gid_ttype$biotype[idx]
        
        ensg$gene_id[mask3] <- str_remove(parent_vec, "^gene:")
        ensg$transcript_type[mask3] <- ttype_vec
        ensg$gene_type[mask3] <- unname(gid_to_gtype_vec[ensg$gene_id[mask3]])
    }

    # Add in additional versions of tid and gid for parse_GFF.R
    ensg$transcript_name <- ensg$transcript_id
    ensg$gene_name <- ensg$gene_id

    # Remove redundant columns.
    ensg <- ensg %>% select(-id, -biotype, -parent)
    
    # Create a separate dataphase with only the CDS regions.
    cds <- ensg %>% filter(type == "CDS")
    
    # Define which columns to keep in the dataphase (used by group_by)
    to_keep <- c('seqid','source','score','strand','gene_id','transcript_id',
               'gene_type','transcript_type','transcript_name','gene_name', 'transcript_support_level')
    
    # Add start/stop codons on both strands
    pos_cds_starts <- add_stop_start(cds, to_keep, strand = "+", start = TRUE)
    pos_cds_stops  <- add_stop_start(cds, to_keep, strand = "+", start = FALSE)
    neg_cds_starts <- add_stop_start(cds, to_keep, strand = "-", start = TRUE)
    neg_cds_stops  <- add_stop_start(cds, to_keep, strand = "-", start = FALSE)

    # Concatenate starts/stops with the original table
    all_concat <- bind_rows(pos_cds_starts, pos_cds_stops, neg_cds_starts, neg_cds_stops, ensg)
    
    # Build accession ranking from gene_type present in all_concat
    for_accession <- tibble(accession = unique(all_concat$gene_type))
    master <- readr::read_tsv(master_ranking_path, show_col_types = FALSE)
    acession_ranking <- left_join(for_accession, master, by = "accession") %>%
        filter(!is.na(rank)) %>%
            arrange(rank) %>%
                mutate(rank = seq_len(n()))
    
    readr::write_tsv(acession_ranking, output_accession)
    
    # Perform the conversion: collapse type_2 to "transcript"
    all_concat$type[all_concat$type %in% type_2] <- "transcript"
    
    # Add the chr prefix
    all_concat$seqid <- paste0("chr", all_concat$seqid)
    
    # Export as GFF3
    write_df_as_gff3(all_concat, output_gff)
}

# Load in the gff file. 
df <- read_gff_as_df(input_gff)

if (tolower(source_mode) == "ensembl") {
  df <- df |> as.data.frame()
  process_ensembl(df, master_ranking_path, output_accession, output_gff)
} else if (tolower(source_mode) == "gencode") {
  process_gencode(df, master_ranking_path, output_accession, output_gff)
} else {
  stop('Unsupported --source: ', source_mode, ' (expected "ensembl" or "gencode")')
}
