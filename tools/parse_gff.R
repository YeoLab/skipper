#!/home/eboyle/bin/Rscript --vanilla 

library(magrittr)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
gff3_file = args[1]
accession_ranking_file = args[2]
partition_output = args[3]
annotations_output = args[4]
window_size = 100

transcript_data = rtracklayer::readGFF(gff3_file) %>% .[which(.$transcript_type != "artifact"),]
gr = makeGRangesFromDataFrame(transcript_data, keep.extra.columns=TRUE)
accession_data = readr::read_tsv(accession_ranking_file) %>% (dplyr::arrange)(rank)
accession_type_rankings = c(accession_data$accession, "primary_miRNA")
exon_subtypes = accession_data$exon_subtype %>% unique
protein_coding_subtype = accession_data$exon_subtype[accession_data$accession == "protein_coding"] %>% head(1)
prioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) < 1]
unprioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) >= 1]

feature_order = c(paste0("EXON_", prioritized_exon_subtypes), "CDS_SOLITARY", "CDS_START","CDS_STOP","CDS","UTR3","UTR5",paste0("EXON_", unprioritized_exon_subtypes),"SSB_ADJ","SS3_ADJ","SS5_ADJ","SSB_PROX","SS3_PROX","SS5_PROX","PRIMIRNA","INTRON")
write("Checking that all gene and transcript types are included in ranking", stderr())
stopifnot(length(setdiff(unique(gr$gene_type), accession_type_rankings)) == 0, length(setdiff(unique(gr$transcript_type), accession_type_rankings)) == 0)
write("...Success", stderr())

gr$metadata = paste0(gr$gene_name, ":", gr$gene_id, ":", gr$transcript_id, ":", gr$gene_type, ":", gr$transcript_type)

# UTRs
utr3 = gr[(gr$type == "three_prime_UTR")] %>% sort
utr5 = gr[(gr$type == "five_prime_UTR")] %>% sort

# CDSs
cds = gr[gr$type == "CDS"] %>% sort
cds_split = split(cds, cds$metadata)

start_codon = gr[gr$type == "start_codon"]
stop_codon = gr[gr$type == "stop_codon"]
cds_start = findOverlaps(cds, start_codon) %>% (function(overlaps) cds[queryHits(overlaps[cds$metadata[queryHits(overlaps)] == start_codon$metadata[subjectHits(overlaps)]])]) %>% resize(.,ifelse(width(.) < 100, width(.), 100)) %>% sort
cds_stop = findOverlaps(cds, stop_codon) %>% (function(overlaps) cds[queryHits(overlaps[cds$metadata[queryHits(overlaps)] == stop_codon$metadata[subjectHits(overlaps)]])]) %>% resize(.,ifelse(width(.) < 100, width(.), 100),fix="end") %>% sort

cds_start_split = split(cds_start, cds_start$metadata)
cds_stop_split = split(cds_stop, cds_stop$metadata)
first_last_transcripts = intersect(names(cds_start_split), names(cds_stop_split))
cds_single = intersect(cds_start_split[first_last_transcripts], cds_stop_split[first_last_transcripts]) %>% stack("metadata") %>% sort 

# miRNA padding
primirna = (gr[which(gr$transcript_type == "miRNA")] + 500) %>% sort
if (length(primirna) > 0){                                    
    primirna$transcript_type = "primary_miRNA" 
    primirna$metadata = paste0(primirna$gene_name, ":", primirna$gene_id, ":", primirna$transcript_id, ":", primirna$gene_type, ":", primirna$transcript_type)
} else {
    primirna = NULL
}

# exons (not CDS)
exons = gr[gr$type == "exon"] %>% sort
exons_split = split(exons, exons$metadata) 

transcripts = gr[gr$type == "transcript"] %>% sort
transcripts_split = split(transcripts, transcripts$metadata)

# introns
introns_split = setdiff(transcripts_split, exons_split)
introns_unsplit = introns_split %>% stack("metadata")

# splice site adjacent
ss5_adj = resize(introns_unsplit,ifelse(width(introns_unsplit) < 100, width(introns_unsplit), 100)) 
ss3_adj = resize(introns_unsplit,ifelse(width(introns_unsplit) < 100, width(introns_unsplit), 100), fix = "end") 
ssb_adj = intersect(split(ss5_adj,ss5_adj$metadata), split(ss3_adj,ss3_adj$metadata)) %>% stack("metadata")

# splice site proximal
ss5_prox = resize(introns_unsplit,ifelse(width(introns_unsplit) < 500, width(introns_unsplit), 500)) 
ss3_prox = resize(introns_unsplit,ifelse(width(introns_unsplit) < 500, width(introns_unsplit), 500), fix = "end") 
ssb_prox = intersect(split(ss5_prox,ss5_prox$metadata), split(ss3_prox,ss3_prox$metadata)) %>% stack("metadata")

reduce_grange <- function(data, f_type) {
  # Return an empty GRanges with the right column when there is no data.
  if (length(data) == 0L) {
    gr0 <- GRanges()
    mcols(gr0)$feature_type <- character(0)
    return(gr0)
  }
  split_reduction <- (tidyr::separate((dplyr::as_tibble)(data), metadata,
                                      c("gene_name","gene_id","transcript_id","gene_type","transcript_type"), ":")) %>%
    GRanges() %>% split(., .$transcript_type) %>% reduce
  feature_partition <- GRanges()
  for (transcript_type in intersect(accession_type_rankings, names(split_reduction))) {
    subtracted_data <- setdiff(split_reduction[[transcript_type]], feature_partition)
    feature_partition <- c(feature_partition, subtracted_data)
  }
  feature_partition$feature_type <- f_type
  feature_partition
}

# Add subtype to exon features
exons$exon_subtype = dplyr::tibble(transcript_type=exons$transcript_type) %>% dplyr::left_join(., dplyr::rename(accession_data, transcript_type = accession)) %>% dplyr::pull(exon_subtype)

exons_reduced_list = lapply(exon_subtypes, function(subtype) reduce_grange(exons[exons$exon_subtype == subtype], paste0("EXON_", subtype)))
prioritized_exon_subtypes_reduced = exons_reduced_list[cumsum(exon_subtypes == protein_coding_subtype) < 1]
unprioritized_exon_subtypes_reduced = exons_reduced_list[cumsum(exon_subtypes == protein_coding_subtype) >= 1]

cds_single_reduced = reduce_grange(cds_single, "CDS_SOLITARY")
cds_start_reduced = reduce_grange(cds_start, "CDS_START") 
cds_stop_reduced = reduce_grange(cds_stop, "CDS_STOP") 
cds_reduced = reduce_grange(cds, "CDS") 
utr3_reduced = reduce_grange(utr3, "UTR3")
utr5_reduced = reduce_grange(utr5, "UTR5")
## exons are typically broken down by subtype
# exons_reduced = reduce_grange(exons, "EXON")
ssb_adj_reduced = reduce_grange(ssb_adj, "SSB_ADJ")
ss3_adj_reduced = reduce_grange(ss3_adj, "SS3_ADJ")
ss5_adj_reduced = reduce_grange(ss5_adj, "SS5_ADJ")
ssb_prox_reduced = reduce_grange(ssb_prox, "SSB_PROX")
ss3_prox_reduced = reduce_grange(ss3_prox, "SS3_PROX")
ss5_prox_reduced = reduce_grange(ss5_prox, "SS5_PROX")                               
primirna_reduced = reduce_grange(primirna, "PRIMIRNA")
introns_reduced = reduce_grange(introns_unsplit, "INTRON")

feature_data_list <- c(
  prioritized_exon_subtypes_reduced,
  list(
    cds_single_reduced,
    cds_start_reduced,
    cds_stop_reduced,
    cds_reduced,
    utr3_reduced,
    utr5_reduced
  ),
  unprioritized_exon_subtypes_reduced,
  list(
    ssb_adj_reduced,
    ss3_adj_reduced,
    ss5_adj_reduced,
    ssb_prox_reduced,
    ss3_prox_reduced,
    ss5_prox_reduced,
    primirna_reduced,   # This will be empty if there are no miRNAs.
    introns_reduced
  )
)

full_partition <- GRanges()
for (feature_data in feature_data_list) {
  subtracted_data <- setdiff(feature_data, full_partition)
  full_partition <- c(full_partition, subtracted_data)
}

# Tile to fixed windows.
tiled_partition <- sort(full_partition) %>% tile(width = window_size) %>% stack("feature_id")

# Feature annotations (types).
feature_data_concat <- do.call("c", feature_data_list)
feature_hits <- findOverlaps(tiled_partition, feature_data_concat)

feature_annotations <- dplyr::tibble(
  row_id = queryHits(feature_hits),
  feature_type = feature_data_concat$feature_type[subjectHits(feature_hits)]
) %>%
  dplyr::group_by(row_id) %>%
  dplyr::summarize(
    feature_types = stringr::str_flatten(intersect(feature_order, feature_type), collapse = ":"),
    feature_type_top = sub(":.*", "", feature_types),
    .groups = "drop"
  )

# Meta annotations: always concatenate transcripts and primirna; empty primirna is harmless.
meta_hits <- findOverlaps(tiled_partition, c(transcripts, primirna))

meta_annotations <- dplyr::tibble(
  row_id = queryHits(meta_hits),
  metadata = c(transcripts, primirna)$metadata[subjectHits(meta_hits)]
) %>%
  (tidyr::separate)(metadata, c("gene_name","gene_id","transcript_id","gene_type","transcript_type"), sep = ":") %>%
  (dplyr::group_by)(row_id) %>%
  (dplyr::summarize)(
    gene_name = stringr::str_flatten(unique(gene_name), collapse = ":"),
    gene_id = stringr::str_flatten(unique(gene_id), collapse = ":"),
    transcript_ids = stringr::str_flatten(unique(transcript_id), collapse = ":"),
    gene_type_top = intersect(accession_type_rankings, gene_type) %>% head(1),
    transcript_type_top = intersect(accession_type_rankings, transcript_type) %>% head(1),
    gene_types = stringr::str_flatten(intersect(accession_type_rankings, unique(gene_type)), collapse = ":"),
    transcript_types = stringr::str_flatten(intersect(accession_type_rankings, unique(transcript_type)), collapse = ":"),
    .groups = "drop"
  )

tiled_partition$name <- seq_len(length(tiled_partition))
rtracklayer::export.bed(tiled_partition, partition_output)

annotated_features <- dplyr::bind_cols(
  meta_annotations,
  feature_annotations %>% dplyr::select(-row_id)
) %>%
  (dplyr::mutate)(feature_id = tiled_partition$feature_id %>% as.numeric()) %>%
  dplyr::group_by(feature_id) %>% dplyr::mutate(feature_bin = dplyr::row_number()) %>% (dplyr::ungroup) %>%
  dplyr::mutate(
    chrom = as.character(seqnames(tiled_partition)),
    start = start(tiled_partition) - 1L,
    end = end(tiled_partition),
    strand = as.character(strand(tiled_partition))
  ) %>%
  dplyr::transmute(
    chrom, start, end, name = row_id, score = 0, strand, feature_id, feature_bin,
    feature_type_top, feature_types, gene_name, gene_id, transcript_ids,
    gene_type_top, transcript_type_top, gene_types, transcript_types
  )

readr::write_tsv(annotated_features, annotations_output)
