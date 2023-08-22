
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

# list all files
gtf_list <-
  list.files(here::here("results"), pattern = "*.gtf", full.names = TRUE)

tsv_list <-
  list.files(here::here("results"), pattern = "*.tsv", full.names = TRUE)

# Since the panel data files might be in the same directory with the same extension they need to be excluded
tsv_list <-
  tsv_list[!grepl(paste(c("gene_set.tsv", "panel_data.tsv"), collapse = "|"), tsv_list)]

# Load data ---------------------------------------------------------------

# Sample information
sample_info <- read.table(here::here("raw_data", "samples.txt"), header = T, sep = "\t")

# Create a list with all the GTF files. GTF will be read a GRanges by rtracklayer
sample_gtf <-
  lapply(gtf_list, function(x) {
    rtracklayer::import(x) 
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

uniqueTranscripts <- function(data, exon_type) {
  
  # Check if the data argument is a GRanges object or a list of GRanges or a data.frame object
  if (!inherits(data, "GRanges") && !is.list(data) && !is.data.frame(data)) {
    
    stop("data must be either a GRanges object, a list of GRanges objects, or a data.frame object")
    
  }
  
  # Check if the exon_type argument is "exon" or "CDS"
  if (!exon_type %in% c("exon", "CDS")) {
    
    stop("exon_type must be either 'exon' or 'CDS'")
    
  }
  
  # Combine all input files into a single GRanges object or data.frame
  if (inherits(data, "GRanges")) {
    
    all_data <- as.data.frame(data)
    
  } else if (is.list(data)) {
    
    all_data <- do.call(rbind.fill, lapply(data, as.data.frame))
    
  } else if (is.data.frame(data)) {
    
    all_data <- data
    
  }
  
  # this data frame contains all transcripts and will be used
  all_data_distinct <- distinct(all_data)
  
  all_data_filtered <- all_data_distinct[all_data_distinct$type == exon_type,] # Filter the GRanges object or data.frame to contain only exon_type
  all_data_grouped <- split(all_data_filtered, all_data_filtered$transcript_id) # This will split by transcript
  
  result_list <- lapply(all_data_grouped, function(transcript) {
    
    transcript_str <- paste(
      paste(unique(transcript$seqnames), collapse = ","),
      paste(unique(transcript$strand), collapse = ","),
      paste(
        paste(transcript$start, transcript$end, sep = ","), 
        collapse = ","
      ),
      sep = ","
    )
    transcript_id <- unique(transcript$transcript_id)
    return(data.frame(transcript_str = transcript_str, transcript_id = transcript_id))
    
  })
  
  # Combine the list of data frames into a single data frame
  string_df <- do.call(rbind, result_list)
  string_df <- string_df[order(string_df$transcript_str), ] # sort to speed up comparisons
  row.names(string_df) <- NULL
  
  # Comparison using the sorted data
  # will update the code as before this it is all base R
  grouped_data <- group_by(string_df, transcript_str)
  summarized_data <- summarize(grouped_data, all_transcripts = paste(transcript_id, collapse = ","))
  summarized_data <- dplyr::mutate(summarized_data, 
                                   unique_id = str_extract(all_transcripts, "^[^,]*"),
                                   matching_id = all_transcripts) %>% 
    separate_rows(matching_id, sep = ",")
    
  # need to make sure the matching is on the column with duplicated ids
  final_df <- all_data_distinct %>% 
    dplyr::left_join(., dplyr::select(summarized_data, matching_id, unique_id, all_transcripts), by = c("transcript_id" = "matching_id")) %>% 
    dplyr::mutate(transcript_id = unique_id) %>% 
    dplyr::select(-c(unique_id, transcript_name, talon_transcript, exon_id, talon_exon, talon_gene)) %>% # remove these columns as otherwise distinct will not work since they have different ids
    distinct()
  
  return(final_df)
}

# Main --------------------------------------------------------------------

# unique transcripts without any fuzzy mathing
unique_transcripts <- uniqueTranscripts(data = sample_gtf, exon_type = "exon")

# Save data ---------------------------------------------------------------

saveRDS(unique_transcripts, file = here::here("results", "unique_transcripts.rds"))
