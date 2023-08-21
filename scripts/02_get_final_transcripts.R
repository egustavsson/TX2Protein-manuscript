

# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(GenomicRanges)
  library(rtracklayer)
})


# Load data ---------------------------------------------------------------

unique_transcripts <- readRDS(here::here("results", "unique_transcripts.rds"))

expression_RPM <- readRDS(here::here("results", "expression_RPM.rds"))

# Main --------------------------------------------------------------------

# Filter to only include transcripts of interest wit novel ORFs from genes of interest
transcripts_to_include <- 
  expression_RPM %>% 
  dplyr::select(annot_gene_id, annot_transcript_id, transcript_novelty, reads) %>% 
  dplyr::filter(annot_transcript_id %in% unique_transcripts$transcript_id, # only include unique transcripts
                transcript_novelty %in% c("NNC", "NIC"), # only include transcrpts with novel ORFs
                reads > 1) %>% # only include transcripts with at least 2 reads
  distinct() %>% 
  .$annot_transcript_id

# Final transcript file
transcripts_novel_ORF <- 
  unique_transcripts %>% 
  dplyr::filter(transcript_id %in% transcripts_to_include)

# Save data ---------------------------------------------------------------

rtracklayer::export(transcripts_novel_ORF, 
                    here::here("results", "transcripts_novel_ORF.gff"), format = "GFF3")
