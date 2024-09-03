
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

# Main --------------------------------------------------------------------

# Only include genes expressed across all samples
genes_to_include <- 
  expression_gene_set %>% 
  dplyr::select(annot_gene_id, sample) %>% 
  group_by(annot_gene_id) %>%
  summarize(n_samples = n_distinct(sample)) %>%
  dplyr::filter(n_samples == 9)

# Filter to only include transcripts of interest wit novel ORFs from genes of interest
transcripts_to_include <- 
  expression_gene_set %>% 
  dplyr::select(annot_gene_id, unique_id, transcript_novelty, reads) %>% 
  dplyr::filter(annot_gene_id %in% genes_to_include$annot_gene_id,
                transcript_novelty %in% c("NNC", "NIC"), # only include transcripts with novel ORFs
                reads > 1) %>% # only include transcripts with at least 2 reads
  .$unique_id %>% 
  unique()

# Final transcript file
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
