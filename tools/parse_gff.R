# Load required packages.
library(magrittr)
library(GenomicRanges)

# Load up arguments.
args = commandArgs(trailingOnly=TRUE)
gff3_file = args[1]              
accession_ranking_file = args[2]  
partition_output = args[3]        
annotations_output = args[4]     
window_size = 100                 


######################### Read annotation data #############################
# Import GFF3, excluding transcripts of type "artifact".
transcript_data = rtracklayer::readGFF(gff3_file) %>% .[which(.$transcript_type != "artifact"),]
gr = makeGRangesFromDataFrame(transcript_data, keep.extra.columns=TRUE)

# Import accession rankings and define type/subtype priorities.
accession_data = readr::read_tsv(accession_ranking_file) %>% dplyr::arrange(rank)
accession_type_rankings = c(accession_data$accession, "primary_miRNA")
exon_subtypes = accession_data$exon_subtype %>% unique
protein_coding_subtype = accession_data$exon_subtype[accession_data$accession == "protein_coding"] %>% head(1)
prioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) < 1]
unprioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) >= 1]

# Define a canonical order of feature types.
feature_order = c(
    paste0("EXON_", prioritized_exon_subtypes),
    "CDS_SOLITARY", "CDS_START", "CDS_STOP", "CDS", "UTR3", "UTR5",
    paste0("EXON_", unprioritized_exon_subtypes),
    "SSB_ADJ", "SS3_ADJ", "SS5_ADJ", "SSB_PROX", "SS3_PROX", "SS5_PROX", "PRIMIRNA", "INTRON"
)

# Check that gene/transcript types in GFF match rankings.
write("Checking that all gene and transcript types are included in ranking.", stderr())

gff_gene_types = sort(unique(stats::na.omit(gr$gene_type)))
gff_transcript_types = sort(unique(stats::na.omit(gr$transcript_type)))
ranking_types = sort(unique(stats::na.omit(accession_type_rankings)))

missing_gene_types = setdiff(gff_gene_types, ranking_types)
missing_transcript_types = setdiff(gff_transcript_types, ranking_types)
unused_ranking_types = setdiff(ranking_types, union(gff_gene_types, gff_transcript_types))

# Always report unused ranking entries.
if (length(unused_ranking_types) > 0) {
    write("Note: The following accession types are present in the ranking file but not used in this GFF:", stderr())
    write(paste0("    - ", unused_ranking_types), stderr())
}

# Only stop if required types are missing.
if (length(missing_gene_types) > 0 || length(missing_transcript_types) > 0) {

    error_lines = c("Accession ranking file does not match the GFF annotation.", "",
                    paste0("Ranking file: ", accession_ranking_file), paste0("GFF file: ", gff3_file), "")

    if (length(missing_gene_types) > 0) {
        error_lines = c(error_lines, "Add the following missing accession type to the ranking file:",
                        paste0("    - ", missing_gene_types),"")
    }
    stop(paste(error_lines, collapse = "\n"), call. = FALSE)
}

write("...Success", stderr())

# Build metadata string for each feature.
gr$metadata = paste0(gr$gene_name, ":", gr$gene_id, ":", gr$transcript_id, ":", gr$gene_type, ":", gr$transcript_type)


######################### Define feature-level GRanges #########################
# UTRs.
utr3 = gr[(gr$type == "three_prime_UTR")] %>% sort
utr5 = gr[(gr$type == "five_prime_UTR")] %>% sort

# CDS regions and codons.
cds = gr[gr$type == "CDS"] %>% sort
cds_split = split(cds, cds$metadata)

start_codon = gr[gr$type == "start_codon"]
stop_codon = gr[gr$type == "stop_codon"]

# Find CDSs overlapping with start/stop codons and resize them to <=100 bp.
cds_start = findOverlaps(cds, start_codon) %>%
  (function(overlaps) cds[queryHits(overlaps[cds$metadata[queryHits(overlaps)] == start_codon$metadata[subjectHits(overlaps)]])]) %>%
  resize(., ifelse(width(.) < 100, width(.), 100)) %>% sort

cds_stop = findOverlaps(cds, stop_codon) %>%
  (function(overlaps) cds[queryHits(overlaps[cds$metadata[queryHits(overlaps)] == stop_codon$metadata[subjectHits(overlaps)]])]) %>%
  resize(., ifelse(width(.) < 100, width(.), 100), fix="end") %>% sort

cds_start_split = split(cds_start, cds_start$metadata)
cds_stop_split = split(cds_stop, cds_stop$metadata)

# Identify transcripts that have both start and stop codons.
first_last_transcripts = intersect(names(cds_start_split), names(cds_stop_split))
cds_single = intersect(cds_start_split[first_last_transcripts],
                       cds_stop_split[first_last_transcripts]) %>% stack("metadata") %>% sort 

# pad with ±500bp to define "primary miRNA" features.
primirna = (gr[which(gr$transcript_type == "miRNA")] + 500) %>% sort
if (length(primirna) > 0){                                    
    primirna$transcript_type = "primary_miRNA" 
    primirna$metadata = paste0(primirna$gene_name, ":", primirna$gene_id, ":", primirna$transcript_id, ":", primirna$gene_type, ":", primirna$transcript_type)
} else {
    primirna = NULL
}

# Exons (excluding CDS).
exons = gr[gr$type == "exon"] %>% sort
exons_split = split(exons, exons$metadata) 

# Transcripts.
transcripts = gr[gr$type == "transcript"] %>% sort
transcripts_split = split(transcripts, transcripts$metadata)

# Introns = transcript span minus exon span.
introns_split = setdiff(transcripts_split, exons_split)
introns_unsplit = introns_split %>% stack("metadata")

# Splice site adjacent windows (~100bp at ends of introns).
ss5_adj = resize(introns_unsplit, ifelse(width(introns_unsplit) < 100, width(introns_unsplit), 100)) 
ss3_adj = resize(introns_unsplit, ifelse(width(introns_unsplit) < 100, width(introns_unsplit), 100), fix="end") 
ssb_adj = intersect(split(ss5_adj,ss5_adj$metadata), split(ss3_adj,ss3_adj$metadata)) %>% stack("metadata")

# Splice site proximal windows (~500bp at ends of introns).
ss5_prox = resize(introns_unsplit, ifelse(width(introns_unsplit) < 500, width(introns_unsplit), 500)) 
ss3_prox = resize(introns_unsplit, ifelse(width(introns_unsplit) < 500, width(introns_unsplit), 500), fix="end") 
ssb_prox = intersect(split(ss5_prox,ss5_prox$metadata), split(ss3_prox,ss3_prox$metadata)) %>% stack("metadata")


######################### reduce overlapping features #########################

# Helper function. 
reduce_grange = function(data, f_type) {
  # Return an empty GRanges with the right column when no data is available
  if (length(data) == 0L) {
    gr0 = GRanges()
    mcols(gr0)$feature_type = character(0)
    return(gr0)
  }
  # Split by transcript type, reduce overlaps, and prioritize by accession ranking
  split_reduction = (tidyr::separate((dplyr::as_tibble)(data), metadata,
                                      c("gene_name","gene_id","transcript_id","gene_type","transcript_type"), ":")) %>%
    GRanges() %>% split(., .$transcript_type) %>% reduce
  feature_partition = GRanges()
  for (transcript_type in intersect(accession_type_rankings, names(split_reduction))) {
    subtracted_data = setdiff(split_reduction[[transcript_type]], feature_partition)
    feature_partition = c(feature_partition, subtracted_data)
  }
  feature_partition$feature_type = f_type
  feature_partition
}

# Reduce exon subtypes separately
exons$exon_subtype = dplyr::tibble(transcript_type=exons$transcript_type) %>%
  dplyr::left_join(., dplyr::rename(accession_data, transcript_type = accession)) %>%
  dplyr::pull(exon_subtype)

exons_reduced_list = lapply(exon_subtypes, function(subtype) 
  reduce_grange(exons[exons$exon_subtype == subtype], paste0("EXON_", subtype)))
prioritized_exon_subtypes_reduced = exons_reduced_list[cumsum(exon_subtypes == protein_coding_subtype) < 1]
unprioritized_exon_subtypes_reduced = exons_reduced_list[cumsum(exon_subtypes == protein_coding_subtype) >= 1]


# Reduce other features
cds_single_reduced = reduce_grange(cds_single, "CDS_SOLITARY")
cds_start_reduced = reduce_grange(cds_start, "CDS_START") 
cds_stop_reduced = reduce_grange(cds_stop, "CDS_STOP") 
cds_reduced = reduce_grange(cds, "CDS") 
utr3_reduced = reduce_grange(utr3, "UTR3")
utr5_reduced = reduce_grange(utr5, "UTR5")
ssb_adj_reduced = reduce_grange(ssb_adj, "SSB_ADJ")
ss3_adj_reduced = reduce_grange(ss3_adj, "SS3_ADJ")
ss5_adj_reduced = reduce_grange(ss5_adj, "SS5_ADJ")
ssb_prox_reduced = reduce_grange(ssb_prox, "SSB_PROX")
ss3_prox_reduced = reduce_grange(ss3_prox, "SS3_PROX")
ss5_prox_reduced = reduce_grange(ss5_prox, "SS5_PROX")                               
primirna_reduced = reduce_grange(primirna, "PRIMIRNA")
introns_reduced = reduce_grange(introns_unsplit, "INTRON")

######################### Combine all reduced feature sets into a partition #########################
feature_data_list = c(
  prioritized_exon_subtypes_reduced,
  list(cds_single_reduced, cds_start_reduced, cds_stop_reduced, cds_reduced, utr3_reduced, utr5_reduced),
  unprioritized_exon_subtypes_reduced,
  list(ssb_adj_reduced, ss3_adj_reduced, ss5_adj_reduced,
       ssb_prox_reduced, ss3_prox_reduced, ss5_prox_reduced,
       primirna_reduced, introns_reduced)
)

full_partition = GRanges()
for (feature_data in feature_data_list) {
  subtracted_data = setdiff(feature_data, full_partition)
  full_partition = c(full_partition, subtracted_data)
}

# Base features represents the broad biological region without contextual annotations like splice-site proximity. 
base_feature_list = c(
    list(cds_reduced, utr3_reduced, utr5_reduced),
    exons_reduced_list,
    list(primirna_reduced, introns_reduced)
)

base_partition = GRanges()
for (base_feature_data in base_feature_list) {
    subtracted_data = setdiff(base_feature_data, base_partition)
    base_partition = c(base_partition, subtracted_data)
}

base_partition = sort(base_partition)
base_partition$base_feature_id = seq_len(length(base_partition))
                            
# Tile reduced features into fixed-size windows
tiled_partition = sort(full_partition) %>% tile(width = window_size) %>% stack("feature_id")

######################### Annotate tiled windows by overlapping features #########################
feature_data_concat = do.call("c", feature_data_list)
feature_hits = findOverlaps(tiled_partition, feature_data_concat)

feature_annotations = dplyr::tibble(
  row_id = queryHits(feature_hits),
  feature_type = feature_data_concat$feature_type[subjectHits(feature_hits)]
) %>%
  dplyr::group_by(row_id) %>%
  dplyr::summarize(
    feature_types = stringr::str_flatten(intersect(feature_order, feature_type), collapse = ":"),
    feature_type_top = sub(":.*", "", feature_types),
    .groups = "drop"
  )

# Annotate tiled windows by base feature.
base_hits = findOverlaps(tiled_partition, base_partition)

base_annotations = dplyr::tibble(
    row_id = queryHits(base_hits),
    base_feature_type = base_partition$feature_type[subjectHits(base_hits)],
    base_feature_id = base_partition$base_feature_id[subjectHits(base_hits)]
) %>%
    dplyr::distinct(row_id, .keep_all = TRUE)

# Meta-annotations (gene/transcript IDs).
meta_features = c(transcripts, primirna)
meta_hits = findOverlaps(tiled_partition, meta_features)

# Convert overlap result to a data.table using metadata columns directly.
meta_dt = data.table::data.table(
    row_id = queryHits(meta_hits),
    gene_name = meta_features$gene_name[subjectHits(meta_hits)],
    gene_id = meta_features$gene_id[subjectHits(meta_hits)],
    transcript_id = meta_features$transcript_id[subjectHits(meta_hits)],
    gene_type = meta_features$gene_type[subjectHits(meta_hits)],
    transcript_type = meta_features$transcript_type[subjectHits(meta_hits)]
)

# Helper to mimic original behavior with data table optimizations.
get_top_ranked_type = function(x, rankings) {
    hits = rankings[rankings %in% x]
    if (length(hits) == 0L) {
        return(NA_character_)
    }
    hits[[1]]
}

# Perform str flatten on everything. 
meta_annotations = meta_dt[
    ,
    .(
        gene_name = stringr::str_flatten(unique(gene_name), collapse = ":"),
        gene_id = stringr::str_flatten(unique(gene_id), collapse = ":"),
        transcript_ids = stringr::str_flatten(unique(transcript_id), collapse = ":"),
        gene_type_top = get_top_ranked_type(gene_type, accession_type_rankings),
        transcript_type_top = get_top_ranked_type(transcript_type, accession_type_rankings),
        gene_types = stringr::str_flatten(
            accession_type_rankings[accession_type_rankings %in% unique(gene_type)],
            collapse = ":"
        ),
        transcript_types = stringr::str_flatten(
            accession_type_rankings[accession_type_rankings %in% unique(transcript_type)],
            collapse = ":"
        )
    ),
    by = row_id
]

# Convert back to tibble/data.frame. 
meta_annotations = tibble::as_tibble(meta_annotations)

######################### Export partition (BED) and annotations (TSV) #########################
tiled_partition$name = seq_len(length(tiled_partition))
rtracklayer::export.bed(tiled_partition, partition_output)

annotated_features = meta_annotations %>%
    dplyr::left_join(feature_annotations, by = "row_id") %>%
    dplyr::left_join(base_annotations, by = "row_id") %>%
    dplyr::mutate(
        feature_id = as.numeric(tiled_partition$feature_id)
    ) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::mutate(feature_bin = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(base_feature_id) %>%
    dplyr::mutate(base_feature_bin = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        chrom = as.character(seqnames(tiled_partition)),
        start = start(tiled_partition) - 1L,
        end = end(tiled_partition),
        strand = as.character(strand(tiled_partition))
    ) %>%
    dplyr::transmute(
        chrom, start, end, name = row_id, score = 0, strand,
        feature_id, feature_bin,
        base_feature_id,
        feature_type_top, feature_types,
        gene_name, gene_id, transcript_ids,
        gene_type_top, transcript_type_top, gene_types, transcript_types
    )

readr::write_tsv(annotated_features, annotations_output)
