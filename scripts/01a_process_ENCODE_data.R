options(warn = -1)

# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(plyr)
  library(here)
  library(GenomicRanges)
  library(plyranges)
})

# Arguments ---------------------------------------------------------------

# List all GTF files
gtf_list <-
  list.files(here::here("results", "ENCODE"), pattern = "*.gtf", full.names = TRUE)

# List all TSV expression data files
tsv_list <-
  list.files(here::here("raw_data"), pattern = "*.tsv", full.names = TRUE)

# Load genes of interest
genes_of_interest <- read.table(here::here("results", "gene_set.tsv"), header = F, sep = "\t", col.names = "genes")

# Load data ---------------------------------------------------------------

# Sample information
sample_info <- read.table(here::here("raw_data", "samples.txt"), header = T, sep = "\t")

# Create a list with all the GTF files. GTF will be read as GRanges by rtracklayer
sample_gtf <-
  lapply(gtf_list, function(x) {
    rtracklayer::import(x) 
  })

# Functions ---------------------------------------------------------------

# Function to extract unique transcripts
uniqueTranscripts <- function(data, exon_type) {
  
  # Check if the data argument is valid
  if (!inherits(data, "GRanges") && !is.list(data) && !is.data.frame(data)) {
    stop("data must be either a GRanges object, a list of GRanges objects, or a data.frame object")
  }
  
  # Check exon_type
  if (!exon_type %in% c("exon", "CDS")) {
    stop("exon_type must be either 'exon' or 'CDS'")
  }
  
  # Combine all input files into a single data frame
  all_data <- if (inherits(data, "GRanges")) {
    as.data.frame(data)
  } else if (is.list(data)) {
    do.call(rbind.fill, lapply(data, as.data.frame))
  } else {
    data
  }
  
  all_data_distinct <- distinct(all_data)
  all_data_filtered <- all_data_distinct[all_data_distinct$type == exon_type, ]
  all_data_grouped <- split(all_data_filtered, all_data_filtered$transcript_id)
  
  result_list <- lapply(all_data_grouped, function(transcript) {
    transcript_str <- paste(
      paste(unique(transcript$seqnames), collapse = ","),
      paste(unique(transcript$strand), collapse = ","),
      paste(paste(transcript$start, transcript$end, sep = ","), collapse = ","),
      sep = ","
    )
    transcript_id <- unique(transcript$transcript_id)
    return(data.frame(transcript_str = transcript_str, transcript_id = transcript_id))
  })
  
  string_df <- do.call(rbind, result_list)
  string_df <- string_df[order(string_df$transcript_str), ]
  row.names(string_df) <- NULL
  
  grouped_data <- group_by(string_df, transcript_str)
  summarized_data <- summarize(grouped_data, all_transcripts = paste(transcript_id, collapse = ","))
  summarized_data <- dplyr::mutate(summarized_data, 
                                   unique_id = str_extract(all_transcripts, "^[^,]*"),
                                   matching_id = all_transcripts) %>% 
    separate_rows(matching_id, sep = ",")
  
  final_df <- all_data_distinct %>% 
    dplyr::left_join(., dplyr::select(summarized_data, matching_id, unique_id, all_transcripts), by = c("transcript_id" = "matching_id")) %>% 
    dplyr::mutate(transcript_id = unique_id) %>% 
    dplyr::select(-c(unique_id, 
                     transcript_name, 
                     talon_transcript, 
                     exon_id, 
                     talon_exon, 
                     talon_gene,
                     ISM.prefix_to_IDs,
                     ISM.prefix_transcript,
                     ISM_to_IDs,
                     ISM.suffix_to_IDs,
                     ISM.suffix_transcript,
                     protein_id,
                     source.1,
                     tag,
                     transcript_support_level,
                     havana_transcript,
                     genomic_transcript,
                     NNC_transcript,
                     ISM_transcript,
                     NIC_transcript)) %>% 
    distinct()
  
  return(final_df)
}

# Function to calculate RPM
getRPM <- function(expression_data, sample_data) {
  
  expression_list <- lapply(expression_data, function(x){
    read.table(x, header = TRUE, sep = "\t") %>% 
      dplyr::mutate(sample = gsub('.tsv', '', basename(x))) %>% 
      rename_if(startsWith(names(.), "rep"), ~"reads")
  })
  
  expression_df_list <- lapply(expression_list, function(df) {
    total_reads <- sum(df$reads)
    per_million_scaling_factor <- total_reads / 1e6
    df$RPM <- df$reads / per_million_scaling_factor
    return(df)
  })
  
  combined_expression_df <- bind_rows(expression_df_list) %>% 
    left_join(., sample_data, by = c("sample" = "sample"))
  
  return(combined_expression_df)
}

# Main --------------------------------------------------------------------

# Extract unique transcripts from the GTF files
unique_transcripts <- uniqueTranscripts(data = sample_gtf, exon_type = "exon")

# Save unique transcripts
saveRDS(unique_transcripts, file = here::here("results", "unique_transcripts.rds"))

# Process expression data -------------------------------------------------

# Get unique transcript IDs
transcript_IDs <- unique_transcripts %>% 
  dplyr::filter(type == "transcript") %>% 
  dplyr::select(transcript_id, all_transcripts) %>% 
  separate_rows(all_transcripts, sep = ",")

# Calculate RPM for expression data
expression_RPM <- getRPM(expression_data = tsv_list, sample_data = sample_info)

# Filter expression data for genes of interest and add unique transcript IDs
expression_gene_set <-
  expression_RPM %>%
  dplyr::filter(str_extract(annot_gene_id, "[^.]+") %in% genes_of_interest$genes) %>% 
  dplyr::mutate(
    unique_id = case_when(
      annot_transcript_id %in% transcript_IDs$all_transcripts ~ {
        matching_ids <- transcript_IDs$transcript_id[match(annot_transcript_id, transcript_IDs$all_transcripts)]
        ifelse(is.na(matching_ids), NA, matching_ids)
      },
      TRUE ~ NA_character_
    )
  ) %>% 
  dplyr::select(-c(transcript_ID, annot_transcript_name, annot_transcript_id))

# Save filtered expression data
saveRDS(expression_gene_set, file = here::here("results", "expression_gene_set.rds"))
