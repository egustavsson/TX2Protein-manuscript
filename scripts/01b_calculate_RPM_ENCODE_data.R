# !!! This script needs to be run on the expression data pre any filtering for genes of interest

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Arguments ---------------------------------------------------------------

# Expression data from ENCODE
tsv_list <-
  list.files(here::here("raw_data"), pattern = "*.tsv", full.names = TRUE)

# Load data ---------------------------------------------------------------

# Sample information
sample_info <- read.table(here::here("raw_data", "samples.txt"), header = T, sep = "\t")

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

expression_RPM <- getRPM(expression_data = tsv_list,
                         sample_data = sample_info)

# Save data ---------------------------------------------------------------

saveRDS(expression_RPM, file = here::here("results", "expression_RPM.rds"))
