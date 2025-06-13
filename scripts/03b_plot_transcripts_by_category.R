# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(scales)

# Load data ---------------------------------------------------------------

expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))

# Functions ---------------------------------------------------------------

TXfilter <- function(data, gene_expression_sample_count, min_read_count, allowed_novelty_types, min_transcript_occurrences) {
  
  selected_genes <- data %>%
    dplyr::select(annot_gene_id, sample) %>%
    group_by(annot_gene_id) %>%
    summarize(sample_count = n_distinct(sample), .groups = "drop") %>%
    dplyr::filter(sample_count >= gene_expression_sample_count) %>%
    dplyr::pull(annot_gene_id)
  
  selected_transcripts <- data %>%
    dplyr::select(unique_id, sample) %>%
    group_by(unique_id) %>%
    summarize(transcript_occurrences = n_distinct(sample), .groups = "drop") %>%
    dplyr::filter(transcript_occurrences >= min_transcript_occurrences) %>%
    dplyr::pull(unique_id)
  
  data %>%
    dplyr::filter(
      annot_gene_id %in% selected_genes,
      transcript_novelty %in% allowed_novelty_types, 
      reads >= min_read_count,
      unique_id %in% selected_transcripts
    )
}

# Main --------------------------------------------------------------------

filtered_data <- TXfilter(
  data = expression_gene_set,
  gene_expression_sample_count = 9, 
  min_transcript_occurrences = 1,
  min_read_count = 2,
  allowed_novelty_types = c("Known", "ISM", "NNC", "NIC")
)

# Get number of transcript per category for each sample
summary_data <- filtered_data %>%
  group_by(sample, transcript_novelty) %>%
  summarise(n_transcripts = n(), .groups = "drop")

# Get absolute number of unique transcripts per category
summary_data_absolute <- filtered_data %>%
  group_by(transcript_novelty) %>%
  summarise(n_unique_transcripts = n_distinct(unique_id), .groups = "drop")

# fill colour to use
fill_colour <- c("Known" = "#4d9221ff",
                 "ISM" = "#74add1ff",
                 "NIC" = "#d53e4fff",
                 "NNC" = "#b2abd2ff")

# Plot transcripts by category
transcripts_per_sample_plot <-
  summary_data %>% 
  ggplot(
    aes(x = transcript_novelty, y = n_transcripts, fill = transcript_novelty)) +
  geom_violin(scale = "width", show.legend = F, trim = T, adjust = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  scale_fill_manual(values = fill_colour) +
  theme_bw() +
  labs(
    title = "Transcripts by Transcript Category",
    x = "Transcript Category",
    y = "Transcripts per Sample"
  ) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.text = element_text(size = 14, face = "bold"))

# Plot absolute number bar plot
transcripts_absolute_plot <-
  summary_data_absolute %>%
  ggplot(aes(x = transcript_novelty, y = n_unique_transcripts, fill = transcript_novelty)) +
  geom_bar(stat = "identity", show.legend = F, color = "black") +
  geom_text(aes(label = n_unique_transcripts), vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = fill_colour) +
  theme_bw() +
  labs(
    title = "Absolute Number of Unique Transcripts by Category",
    x = "Transcript Category",
    y = "Number of Unique Transcripts"
  ) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.text = element_text(size = 14, face = "bold"))

# Save data ---------------------------------------------------------------

ggsave(
  plot = transcripts_per_sample_plot, 
  filename = "SFig2_transcripts_per_sample.svg", 
  path = here::here("results", "plots"), 
  width = 9, 
  height = 7, 
  dpi = 600, 
  bg = "white"
)

ggsave(
  plot = transcripts_absolute_plot, 
  filename = "03b_sum_transcripts.svg", 
  path = here::here("results", "plots"), 
  width = 10, 
  height = 7, 
  dpi = 600, 
  bg = "white"
)
