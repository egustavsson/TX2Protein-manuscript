# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggtranscript)
  library(readxl)
})

# Arguments ---------------------------------------------------------------

args <- list(
  path_to_TX2P = here::here("results", "TX2P_output", "EPI_Revised_main.xlsx"),
  path_to_expression = here::here("results", "expression_gene_set.rds"),
  path_to_GFF = here::here("results", "transcripts_novel_ORF.gff")
)


# Load data ---------------------------------------------------------------

# TX2P output data
# This might need updating depending on the final output format we decide on for TX2P
# TX2P_count contains the number of times each transcript had a positive peptide detection
# TX2P_merged has all the data
TX2P_out <- list(
  TX2P_count = read_excel(args$path_to_TX2P, sheet = "Count_Found_in_mass_Datasets", col_names = c("data_set", "Transcript_ID", "occurences")),
  TX2P_merged = read_excel(args$path_to_TX2P, sheet = "merged_data"),
  TX2P_ORF = read_excel(args$path_to_TX2P, sheet = "ORFS")
)

# This data is generated through 01b_process_ENCODE_expression.R
expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))

GFF <- rtracklayer::import(args$path_to_GFF)


# Functions ---------------------------------------------------------------

TXfilter <- function(data, gene_expression_sample_count, min_read_count, allowed_novelty_types, min_transcript_occurrences) {
  
  # Step 1: Identify genes expressed in at least gene_expression_sample_count samples
  selected_genes <- 
    data %>%
    dplyr::select(annot_gene_id, sample) %>%
    group_by(annot_gene_id) %>%
    summarize(sample_count = n_distinct(sample), .groups = "drop") %>%
    dplyr::filter(sample_count >= gene_expression_sample_count) %>%
    dplyr::pull(annot_gene_id)
  
  # Step 2: Identify transcripts found in at least min_transcript_occurrences samples
  selected_transcripts <- 
    data %>%
    dplyr::select(unique_id, sample) %>%
    group_by(unique_id) %>%
    summarize(transcript_occurrences = n_distinct(sample), .groups = "drop") %>%
    dplyr::filter(transcript_occurrences >= min_transcript_occurrences) %>%
    dplyr::pull(unique_id)
  
  # Step 3: Filter transcripts of interest
  filtered_transcripts <- 
    data %>%
    dplyr::filter(
      annot_gene_id %in% selected_genes,
      transcript_novelty %in% allowed_novelty_types, 
      reads >= min_read_count,
      unique_id %in% selected_transcripts
    )
  
  return(filtered_transcripts)
}

# Main --------------------------------------------------------------------

filtered_data <- TXfilter(
  data = expression_gene_set,
  gene_expression_sample_count = 9, 
  min_transcript_occurrences = 1,
  min_read_count = 2,
  allowed_novelty_types = c("NNC", "NIC")
)


# get filtered data and annotate with MS support
filtered_data_w_MS_support <- filtered_data %>% 
  dplyr::left_join(
    TX2P_out$TX2P_count %>% 
      dplyr::select(Transcript_ID, occurences) %>% 
      distinct(), 
    by = c("unique_id" = "Transcript_ID")
  ) %>% 
  dplyr::mutate(occurences = tidyr::replace_na(occurences, 0))
