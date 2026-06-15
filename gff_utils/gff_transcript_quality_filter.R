# Load core libraries for data manipulation/IO
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# Parse command-line arguments in the expected order.
args = commandArgs(trailingOnly = TRUE)
source_mode = args[1]
input_gff = args[2]
output_gff = args[3]

# Read GFF3 into a data.frame.
read_gff_as_df = function(path) {
    df = rtracklayer::readGFF(path)
    # Make sure we are working with a plain data.frame.
    df = as.data.frame(df, stringsAsFactors = FALSE)
    return(df)
}

# Create helper function for adding in stop and start codon information (needed for ensembl). 
add_stop_start = function(cds, to_keep, strand = "+", start = TRUE) {
  cds_strand = cds %>% filter(strand == strand)
  if (nrow(cds_strand) == 0) return(cds_strand[0, ])

  if (start && strand == "+") {
    df = cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(start = min(start), .groups = "drop") %>%
      mutate(end = start + 2L, type = "start_codon")
  } else if (!start && strand == "+") {
    df = cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(end = max(end), .groups = "drop") %>%
      mutate(start = end - 2L, type = "stop_codon")
  } else if (!start && strand == "-") {
    df = cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(start = min(start), .groups = "drop") %>%
      mutate(end = start + 2L, type = "stop_codon")
  } else if (start && strand == "-") {
    df = cds_strand %>%
      group_by(across(all_of(to_keep))) %>%
      summarise(end = max(end), .groups = "drop") %>%
      mutate(start = end - 2L, type = "start_codon")
  }

  df = df %>%
    mutate(phase = 0L) %>%
    select(any_of(colnames(cds_strand)))
  return(df)
}

# Convert back to GRanges and export as GFF3.
write_df_as_gff3 = function(df, path) {

    # Drop row names.
    rownames(df) = NULL

    # Remove filtering columns if they somehow made it through.
    helper_cols = c("tag_str", "is_mane", "is_canonical", "is_primary", "is_appris", "is_ccds",
                    "is_basic", "priority", "tie_break_score", "min_priority", "max_tie_break_score")
    df = df[, setdiff(colnames(df), helper_cols), drop = FALSE]

    # Coerce list columns to semicolon-delimited character strings.
    df[] = lapply(df, function(col) {
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

    # ensure proper column types. 
    df$start = as.integer(df$start)
    df$end = as.integer(df$end)
    df$strand = as.character(df$strand)

    # Convert to a GRanges object. 
    gr = makeGRangesFromDataFrame(df, seqnames.field = "seqid", start.field = "start", end.field = "end", 
                                  strand.field = "strand", keep.extra.columns = TRUE, ignore.strand = FALSE)

    # Save and gzip the filtered GFF file. 
    con = gzfile(path, "w")
    rtracklayer::export(gr, con, format = "gff3")
    close(con)
}

# process_gencode
process_gencode = function(df, gff_output) {

    # Extract transcript support levels.
    tsl = as.character(df$transcript_support_level)

    # Some transcript support levels are literally the string "NA".
    is_literal_NA = !is.na(tsl) & tsl == "NA"

    # Convert tsl to numeric where possible.
    tsl_num = rep(NA_real_, length(tsl))
    suppressWarnings({
        tsl_num[!is_literal_NA] = as.numeric(tsl[!is_literal_NA])
    })

    # Keep transcripts with TSL <= 3, plus NAs.
    keep = is_literal_NA | is.na(tsl_num) | (tsl_num <= 3)
    df_filtered = df[keep, , drop = FALSE]
    rownames(df_filtered) = NULL

    # Build transcript-level priority table.
    df_tx = df_filtered %>%
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

            priority = case_when(is_mane ~ 1, is_canonical ~ 2, is_primary ~ 3, is_appris ~ 4,
                                 is_ccds ~ 5, is_basic ~ 6, TRUE ~ 7),

            tie_break_score = case_when(priority == 1 ~ 0, 
                                        priority == 2 ~ as.integer(is_primary) + as.integer(is_appris) + as.integer(is_ccds) + as.integer(is_basic),
                                        priority == 3 ~ as.integer(is_appris) + as.integer(is_ccds) + as.integer(is_basic),
                                        priority == 4 ~ as.integer(is_ccds) + as.integer(is_basic),
                                        priority == 5 ~ as.integer(is_basic), TRUE ~ 0)
        ) %>%
        group_by(gene_id) %>%
        mutate(min_priority = min(priority, na.rm = TRUE)) %>%
        ungroup()

    # Keep all top transcripts after tie-breaking.
    best_tx = df_tx %>%
        filter(priority == min_priority) %>%
        group_by(gene_id) %>%
        mutate(max_tie_break_score = max(tie_break_score, na.rm = TRUE)) %>%
        filter(tie_break_score == max_tie_break_score) %>%
        ungroup()

    # Get selected transcript IDs.
    selected_tx_ids = unique(best_tx$transcript_id)

    # Keep all rows belonging to selected transcripts, plus gene rows for selected genes.
    selected_gene_ids = unique(best_tx$gene_id)

    # Perform the final tag filtration. 
    df_out = df_filtered %>%
        filter(gene_id %in% selected_gene_ids &(type == "gene" | transcript_id %in% selected_tx_ids | (type == "transcript" & transcript_id %in% selected_tx_ids)))

    # Write filtered GFF.
    write_df_as_gff3(df_out, gff_output)
}

process_ensembl = function(df, gff_output){
    colnames(df) = tolower(colnames(df))

    # Define 3 annotation type groups, progressively more fine-grained
    type_1 = c("gene", "ncRNA_gene", "pseudogene")
    type_2 = c("mRNA", "lnc_RNA", "ncRNA", "miRNA", "snRNA", "snoRNA", "rRNA", "pseudogenic_transcript", "V_gene_segment",
                "scRNA", "C_gene_segment", "Y_RNA", "J_gene_segment", "D_gene_segment", "transcript", "unconfirmed_transcript")
    type_3 = c("exon", "CDS", "five_prime_UTR", "three_prime_UTR")
    
    # Build mappings for gene-id to gene-type. 
    gid_to_gtype = df %>%
        select(gene_id, biotype) %>%
            filter(!is.na(gene_id), !is.na(biotype)) %>%
                distinct()
    gid_to_gtype_vec = setNames(gid_to_gtype$biotype, gid_to_gtype$gene_id)
    
    # Build mappings for transcript ID to parent gene id and transcript type. 
    tid_to_gid_ttype = df %>%
        select(transcript_id, parent, biotype) %>%
            filter(!is.na(transcript_id), !is.na(parent), !is.na(biotype)) %>%
                distinct()

    # Subset to only the important columns. 
    ensg = df %>% select(seqid, source, type, start, end, score, strand, phase,
                         id, biotype, gene_id, parent, transcript_id, tag)

    # Pre-allocate columns.
    ensg = ensg %>%
        mutate(gene_type = "NA", transcript_type = "NA", gene_id = as.character(gene_id), transcript_id = as.character(transcript_id))
    
    # Gene-level rows.
    mask1 = ensg$type %in% type_1
    ensg$gene_type[mask1] = ensg$biotype[mask1]
    
    # Transcript rows with parent=gene:*
    mask2 = ensg$type %in% type_2
    if (any(mask2)) {
        gene_ids = str_remove(ensg$parent[mask2], "^gene:")
        ensg$gene_id[mask2] = gene_ids
        ensg$gene_type[mask2] = unname(gid_to_gtype_vec[gene_ids])
        ensg$transcript_type[mask2] = ensg$biotype[mask2]
    }

    # type rows with parent=transcript:*
    mask3 = ensg$type %in% type_3
    if (any(mask3)) {
        transcript_ids = str_remove(ensg$parent[mask3], "^transcript:")
        ensg$transcript_id[mask3] = transcript_ids
        
        # map transcript_id -> (parent, biotype)
        idx = match(transcript_ids, tid_to_gid_ttype$transcript_id)
        parent_vec = tid_to_gid_ttype$parent[idx]
        ttype_vec  = tid_to_gid_ttype$biotype[idx]
        
        ensg$gene_id[mask3] = str_remove(parent_vec, "^gene:")
        ensg$transcript_type[mask3] = ttype_vec
        ensg$gene_type[mask3] = unname(gid_to_gtype_vec[ensg$gene_id[mask3]])
    }
    
    # Add in additional versions of tid and gid for parse_GFF.R
    ensg$transcript_name = ensg$transcript_id
    ensg$gene_name = ensg$gene_id
    
    # Remove redundant columns.
    ensg = ensg %>% select(-id, -biotype, -parent)
    
    # Create a separate dataphase with only the CDS regions.
    cds = ensg %>% filter(type == "CDS")

    # Define which columns to keep in the dataphase (used by group_by)
    to_keep = c('seqid','source','score','strand','gene_id','transcript_id',
               'gene_type','transcript_type','transcript_name','gene_name')
    
    # Add start/stop codons on both strands
    pos_cds_starts = add_stop_start(cds, to_keep, strand = "+", start = TRUE)
    pos_cds_stops  = add_stop_start(cds, to_keep, strand = "+", start = FALSE)
    neg_cds_starts = add_stop_start(cds, to_keep, strand = "-", start = TRUE)
    neg_cds_stops  = add_stop_start(cds, to_keep, strand = "-", start = FALSE)

    # Small function that adds "tag" into start/stop dataframes. 
    # (tag column can be list or character, so requires a seperate function). 
    fix_tag_col = function(df, tag_template) {
        template_is_list = is.list(tag_template)
        if (template_is_list) {
            df$tag = I(rep(list("NA"), nrow(df)))
        } else {
            df$tag = rep("NA", nrow(df))
        }
        return(df)
    }

    # Add the appropriate tag columns. 
    pos_cds_starts = fix_tag_col(pos_cds_starts, ensg$tag)
    pos_cds_stops  = fix_tag_col(pos_cds_stops, ensg$tag)
    neg_cds_starts = fix_tag_col(neg_cds_starts, ensg$tag)
    neg_cds_stops  = fix_tag_col(neg_cds_stops, ensg$tag)

    # Add all of the cds start and stop codons to full dataset. 
    all_concat = bind_rows(pos_cds_starts, pos_cds_stops, neg_cds_starts, neg_cds_stops, ensg)
    
    # Collapse type_2 to "transcript"
    all_concat$type[all_concat$type %in% type_2] = "transcript"

    # Collapse type_1 to "gene"
    all_concat$type[all_concat$type %in% type_1] = "gene"
    
    # Add the chr prefix
    all_concat$seqid = paste0("chr", all_concat$seqid)
    
    # Build transcript-level priority table.
    df_tx = all_concat %>%
        filter(type == "transcript") %>%
        mutate(
            tag_str = vapply(tag, function(x) {
                if (length(x) == 0 || all(is.na(x))) {
                    ""
                } else {
                    paste(as.character(x), collapse = ";")
                }
            }, character(1)),
    
            is_mane = str_detect(tag_str, regex("(^|;)MANE_Select($|;)", ignore_case = TRUE)),
            is_mane_clinical = str_detect(tag_str, regex("(^|;)MANE_Plus_Clinical($|;)", ignore_case = TRUE)),
            is_canonical = str_detect(tag_str, regex("(^|;)ensembl_canonical($|;)", ignore_case = TRUE)),
            is_canon_extended = str_detect(tag_str, regex("(^|;)ens_canon_extended($|;)", ignore_case = TRUE)),
            is_appris = str_detect(tag_str, regex("(^|;)appris_principal_[1-5]($|;)", ignore_case = TRUE)),
            is_ccds = str_detect(tag_str, regex("(^|;)CCDS($|;)", ignore_case = TRUE)),
            is_basic = str_detect(tag_str, regex("(^|;)(basic|gencode_basic)($|;)", ignore_case = TRUE)),
            
            priority = case_when(is_mane ~ 1, is_mane_clinical ~ 2, is_canonical ~ 3, is_appris ~ 4,
                                 is_ccds ~ 5, is_canon_extended ~ 6, is_basic ~ 7,TRUE ~ 8),
            
            tie_break_score = case_when(priority == 1 ~ as.integer(is_canonical) + as.integer(is_appris) + as.integer(is_ccds) + as.integer(is_basic),
                                        priority == 2 ~ as.integer(is_canonical) + as.integer(is_appris) + as.integer(is_ccds) + as.integer(is_basic),
                                        priority == 3 ~ as.integer(is_appris) + as.integer(is_ccds) + as.integer(is_basic),
                                        priority == 4 ~ as.integer(is_ccds) + as.integer(is_basic),
                                        priority == 5 ~ as.integer(is_basic),
                                        priority == 6 ~ as.integer(is_basic), TRUE ~ 0)
        ) %>%
        group_by(gene_id) %>%
        mutate(min_priority = min(priority, na.rm = TRUE)) %>%
        ungroup()
    
    # Keep all top transcripts after tie-breaking.
    best_tx = df_tx %>%
        filter(priority == min_priority) %>%
        group_by(gene_id) %>%
        mutate(max_tie_break_score = max(tie_break_score, na.rm = TRUE)) %>%
        filter(tie_break_score == max_tie_break_score) %>%
        ungroup()
    
    # Get selected transcript IDs.
    selected_tx_ids = unique(best_tx$transcript_id)
    
    # Keep all rows belonging to selected transcripts, plus gene rows for selected genes.
    selected_gene_ids = unique(best_tx$gene_id)
    
    # Perform the final tag filtration. 
    df_out = all_concat %>%
        filter(gene_id %in% selected_gene_ids & (type == "gene" | transcript_id %in% selected_tx_ids | (type == "transcript" & transcript_id %in% selected_tx_ids)))
    
    # Write filtered GFF.
    write_df_as_gff3(df_out, gff_output)
}

# Load in the gff file. 
df = read_gff_as_df(input_gff)

if (tolower(source_mode) == "ensembl") {
  process_ensembl(df, output_gff)
} else if (tolower(source_mode) == "gencode") {
  process_gencode(df, output_gff)
} else {
  stop('Unsupported --source: ', source_mode, ' (expected "ensembl" or "gencode")')
}