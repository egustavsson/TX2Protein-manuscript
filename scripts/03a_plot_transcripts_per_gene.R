# Make a histogram as a stacked bar chart with idfferent transcript categories


# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# Load data ---------------------------------------------------------------

unique_transcripts <- readRDS(here::here("results", "unique_transcripts.rds"))

expression_RPM <- readRDS(here::here("results", "expression_RPM.rds"))

# Main --------------------------------------------------------------------

# Get unique IDs as created by 01a_process_ENCODE_transcripts.R and 02....R
transcript_IDs <- 
  unique_transcripts %>% 
  dplyr::filter(type == "transcript") %>% 
  dplyr::select(transcript_id, all_transcripts)


expression_RPM %>% 
  dplyr::mutate(annot_transcript_id = str_extract(annot_transcript_id, "[^.]+"),
                unique_id = case_when(annot_transcript_id %in% transcript_IDs$all_transcripts)) %>% # to finish
  dplyr::select(-c(gene_ID, transcript_ID, annot_transcript_name))



# Save data ---------------------------------------------------------------

rtracklayer::export(transcripts_novel_ORF, 
                    here::here("results", "transcripts_novel_ORF.gff"), format = "GFF3")
