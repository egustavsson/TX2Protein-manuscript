# !!! This script needs to be run on the expression data pre any filtering for genes of interest

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Arguments ---------------------------------------------------------------

# Expression data from ENCODE
tsv_list <-
  list.files(here::here("raw_data"), pattern = "*.tsv", full.names = TRUE)

# After the RPM has been calculated we want to filter to the genes of interest
genes_of_interest <- read.table(here::here("results", "gene_set.tsv"), header = F, sep = "\t", col.names = "genes")

# Load data ---------------------------------------------------------------

# Sample information
sample_info <- read.table(here::here("raw_data", "samples.txt"), header = T, sep = "\t")

# Unique transcripts from 01a_process_ENCODE_transcripts.R
# Need to this data to convert transcript IDs for the expression data
unique_transcripts <- readRDS(here::here("results", "unique_transcripts.rds"))

# Functions ---------------------------------------------------------------

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

# Get unique IDs as created by 01a_process_ENCODE_transcripts.R 
transcript_IDs <- 
  unique_transcripts %>% 
  dplyr::filter(type == "transcript") %>% 
  dplyr::select(transcript_id, all_transcripts) %>% 
  separate_rows(all_transcripts, sep = ",")

# get expression matrix and add unique IDs
expression_RPM <- 
  getRPM(expression_data = tsv_list,
         sample_data = sample_info)
  
  
# Filter to only include genes of interest and add unique IDs
expression_gene_set <-
  expression_RPM %>%
  dplyr::filter(str_extract(annot_gene_id, "[^.]+") %in% genes_of_interest$genes) %>% 
  dplyr::mutate(
    unique_id = case_when(
      annot_transcript_id %in% transcript_IDs$all_transcripts ~ {
        matching_ids <- transcript_IDs$transcript_id[match(annot_transcript_id, transcript_IDs$all_transcripts)]
        ifelse(is.na(matching_ids), NA, matching_ids)
      },
      TRUE ~ NA_character_  # Default value when no conditions are met
    )
  ) %>% 
  dplyr::select(-c(transcript_ID, annot_transcript_name, annot_transcript_id)) # remove these columns as otherwise distinct will not work since they have different ids

# Save data ---------------------------------------------------------------

saveRDS(expression_gene_set, file = here::here("results", "expression_gene_set.rds"))
