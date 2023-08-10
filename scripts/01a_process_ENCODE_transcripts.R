
## THINGS TO CHANGE ##

# Make this a command line script
# compare_transcript function needs to be ran for each gene individually (perhaps a loop or lapply)
# add a step to filter transcripts based on reads and number of samples it is found in
# add an additional script that plots the outcome of this



# Load libraries ----------------------------------------------------------

library(tidyverse)
library(plyr)
library(here)
library(GenomicRanges)
library(plyranges)

# Arguments ---------------------------------------------------------------

# list all files
gtf_list <- 
  list.files(here::here("raw_data"), pattern = "*.gtf", full.names = TRUE)

tsv_list <- 
  list.files(here::here("raw_data"), pattern = "*.tsv", full.names = TRUE)

# Load data ---------------------------------------------------------------

# Sample information
sample_info <- read.table(here::here("raw_data", "samples.txt"), header = T, sep = "\t")

# Create a list with all the GTF files. GTF will be read a GRanges by rtracklayer
sample_gtf <- 
  lapply(gtf_list, function(x) {
    rtracklayer::import(x) %>% 
      plyranges::mutate(sample = gsub('.gtf', '', basename(x)))
  }
  )

# Expression data
expression_data <- 
  lapply(tsv_list, function(x) {
    
    read.table(x, header = TRUE, sep = "\t") %>% 
      dplyr::mutate(sample = gsub('.tsv', '', basename(x))) %>% 
      rename_if(startsWith(names(.), "rep"), ~"reads")
    
  }) %>%
  
  bind_rows(!!!.) %>% 
  left_join(., sample_info, by = c("sample" = "sample")) 

# Functions ---------------------------------------------------------------

# This function compares transcript structures across samples.
compare_transcripts <- function(data, gene_id, exon_type, start_difference = 0, end_difference = 0) {
  
  # Check if the data argument is a GRanges object or a list of GRanges or a data.frame object
  if (!inherits(data, "GRanges") && !is.list(data) && !is.data.frame(data)) {
    
    stop("data must be either a GRanges object, a list of GRanges objects, or a data.frame object")
    
  }
  
  # Check if the gene_id argument is a character vector
  if (!is.character(gene_id)) {
    
    stop("gene_id must be a character vector")
    
  }
  
  # Check if the exon_type argument is "exon" or "CDS"
  if (!exon_type %in% c("exon", "CDS")) {
    
    stop("exon_type must be either 'exon' or 'CDS'")
    
  }
  
  # Check if the start_difference and end_difference arguments are integers
  start_difference <- as.integer(start_difference)
  if (is.na(start_difference)) {
    
    stop("start_difference must be an integer")
    
  }
  end_difference <- as.integer(end_difference)
  if (is.na(end_difference)) {
    
    stop("end_difference must be an integer")
    
  }
  
  # Combine all input files into a single GRanges object or data.frame
  if (inherits(data, "GRanges")) {
    
    all_data <- data
    
  } else if (is.list(data)) {
    
    all_data <- GRanges(do.call(rbind.fill, lapply(data, as.data.frame)))
    
  } else if (is.data.frame(data)) {
    
    all_data <- GRanges(data)
    
  }
  
  # Filter the GRanges object or data.frame to contain only exon_type and specified gene IDs
  all_data_filtered <- all_data[all_data$type == exon_type & all_data$gene_id %in% gene_id,]
  all_data_grouped <- split(all_data_filtered, all_data_filtered$transcript_id)
  
  # Initialize a list to store matching transcript IDs
  matches <- list()
  
  # Precompute values to speed up for loops
  all_names <- names(all_data_grouped)
  all_lengths <- sapply(all_data_grouped, length)
  all_starts <- lapply(all_data_grouped, start)
  all_ends <- lapply(all_data_grouped, end)
  
  # Loop through each transcript
  for (i in all_names) {
    matches[[i]] <- c(i) # include the current transcript as a match
    
    # Loop through each transcript again
    for (j in seq_along(all_names)) {
      
      # if (i == j) next # skip if i and j are the same
      
      if (all_lengths[i] == all_lengths[j]) {
        if (all(rle(countOverlaps(all_data_grouped[[i]], all_data_grouped[[j]], type = "equal"))$values == 1)) {
          if (all(abs(all_starts[[i]] - all_starts[[j]]) <= abs(start_difference)) &
              all(abs(all_ends[[i]] - all_ends[[j]]) <= abs(end_difference))) {
            matches[[i]] <- c(matches[[i]], all_names[j])
          }
        }
      }
    }
  }
  
  # remove duplicated value within each vector
  matches_clean <- lapply(matches, unique)
  
  return(matches_clean[!duplicated(lapply(matches_clean, function(x) paste(sort(x), collapse = "-")))])
  
}

# Main --------------------------------------------------------------------

# Get all unique transcripts and their matches. I have not allowed 5' or 3' differences here
unique_transcripts <- compare_transcripts(data = sample_gtf, 
                                          gene_id = gene_of_interest, exon_type = "exon", 
                                          start_difference = 0, 
                                          end_difference = 0)

# lookup index so we can keep unique transcripts and aggregate
lookup <- stack(unique_transcripts) %>% setNames(c("original_id", "unique_id"))

# get reads per transcript so that we can filter by min coverage. This only includes unique ids
transcript_expression <-
  expression_data %>% 
  left_join(., lookup, by = c("annot_transcript_id" = "original_id")) %>% 
  aggregate(data = ., reads ~ sample + status + transcript_novelty + unique_id, FUN = sum) %>% 
  dplyr::filter(!transcript_novelty %in% c("Genomic"),
                reads >= 2) %>% # here we filter to only include transcripts with at least 2 reads
  group_by(sample) %>% 
  dplyr::mutate(norm_reads = reads/sum(reads) * 100) %>% 
  ungroup()

# Filter gtf to onlyinclude non-matching transcripts
unique_transcript_gtf <-
  sample_gtf %>% 
  GRangesList() %>% # cannot unlist a list of granges to a granges object and therefore need to convert to grangeslist first
  unlist() %>% 
  plyranges::filter(transcript_id %in% unique(transcript_expression$unique_id)) %>% 
  unique() # unique as transcripts are found in  multiple samples and therefore found in multiple

# Save data ---------------------------------------------------------------

saveRDS(transcript_expression, file = here::here("results", "transcript_expression.rds"))
saveRDS(unique_transcript_gtf, file = here::here("results", "unique_transcript_gtf.rds"))

