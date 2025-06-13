# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(GenomicRanges)
  library(rtracklayer)
  library(R.utils)
})


# Functions ---------------------------------------------------------------

# Function to download GENCODE references and import as GRanges

download_and_process_reference <- function(version, genes_of_interest) {
  # Define the file paths dynamically based on the version
  ref_path <- here::here(tempdir(), paste0("gencode.v", version, ".annotation.gtf.gz"))
  unzipped_ref_path <- stringr::str_remove(ref_path, "\\.gz")
  
  # Check if the compressed file exists; if not, download it
  if (!file.exists(ref_path)) {
    download.file(
      url = paste0(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_", version, "/",
        "gencode.v", version, ".annotation.gtf.gz"
      ),
      destfile = ref_path
    )
  }
  
  # Decompress the file
  R.utils::gunzip(ref_path, remove = TRUE)
  
  # Import the GTF file
  ref <- rtracklayer::import(unzipped_ref_path)
  
  # Extract gene IDs from genes_of_interest, ignoring the part after the dot
  gene_ids <- unique(sub("\\..*", "", genes_of_interest$genes))
  
  # Filter reference
  ref_filtered <- ref %>% subset(!type %in% c("gene", "transcript") & sub("\\..*", "", gene_id) %in% gene_ids)
  
  # Remove the uncompressed file to clean up the temporary directory
  file.remove(unzipped_ref_path)
  
  # Return the imported object
  return(ref_filtered)
}


# Function to get unique ranges of a query GRanges compared to a reference GTF

get_unique_ranges <- function(query_gtf, ref_gtf) {
  
  # Reduce and filter query_gtf in one step
  reduced_query_gtf <- query_gtf %>%
    subset(!type %in% c("gene", "transcript")) %>%
    GenomicRanges::reduce()
  
  # Filter ref_gtf for "gene" and "transcript" types
  filtered_ref_gtf <- ref_gtf %>%
    subset(!type %in% c("gene", "transcript")) %>% 
    GenomicRanges::reduce()
  
  # Find overlapping ranges between reduced_query_gtf and filtered_ref_gtf
  overlaps <- findOverlaps(reduced_query_gtf, filtered_ref_gtf)
  
  # Get the ranges from reduced_query_gtf that overlap with filtered_ref_gtf
  overlapping_ranges <- pintersect(
    reduced_query_gtf[queryHits(overlaps)], 
    filtered_ref_gtf[subjectHits(overlaps)]
  )
  
  # Subtract overlapping portions from reduced_query_gtf
  unique_ranges <- GenomicRanges::setdiff(reduced_query_gtf, overlapping_ranges)
  
  # Return the unique ranges
  return(unique_ranges)
}


# Load data ---------------------------------------------------------------

# These are all the epilepsy genes
genes_of_interest <- read.table(here::here("results", "gene_set.tsv"), header = F, sep = "\t", col.names = "genes")

# These are all the transcripts with MS support. This was passed on by David (see email)
final_transcripts <- read_tsv(here::here("results", "TX2P_output", "final_list_of_transcripts.txt"), col_names = F)

# Main --------------------------------------------------------------------

# Get GENCODE references
ref_v29 <- download_and_process_reference(version = 29, genes_of_interest = genes_of_interest)


# Custom transcriptome annotation (e.g., epilepsy annotation)
unique_transcripts <- 
  readRDS(here::here("results", "unique_transcripts.rds")) %>% 
  dplyr::filter(
    transcript_id %in% final_transcripts$X1) %>% 
  GRanges()



# The reference here should be the epilepsy annotation
ref_diff <- 
  get_unique_ranges(
    query_gtf = unique_transcripts, 
    ref_gtf = ref_v29)


sum(width(ref_diff))

sum(width(ref_diff)) / sum(width(ref_v29)) * 100




