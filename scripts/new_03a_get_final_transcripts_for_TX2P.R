
# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(GenomicRanges)
  library(rtracklayer)
})


# Load data ---------------------------------------------------------------

unique_transcripts <- readRDS(here::here("results", "unique_transcripts.rds"))

expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))

# This is after annotating transcpts using pigeon. There seems to be some Talon annotations not being accurate
pigeon <- df <- read.table(here::here("results", "pigeon", "final.tsv"), header = TRUE, sep = "\t")[, c(1, 6)]

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

# Transcript classification done by Pigeon to confirm novelty.
expression_gene_set %>% dplyr::left_join(., pigeon, by = c("unique_id" = "isoform"))





transcripts_to_include <- TXfilter(
  data = expression_gene_set,
  gene_expression_sample_count = 9, 
  min_transcript_occurrences = 1,
  min_read_count = 2,
  allowed_novelty_types = c("NNC", "NIC")
)

# Final transcript GFF file
transcripts_novel_ORF <- 
  unique_transcripts %>% 
  dplyr::filter(transcript_id %in% transcripts_to_include)

# Save data ---------------------------------------------------------------

rtracklayer::export(transcripts_novel_ORF, 
                    here::here("results", "transcripts_novel_ORF.gff"), format = "GFF3")

write.table(data.frame(transcripts_to_include), 
            file = here::here("results", "transcripts.tsv"), sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
