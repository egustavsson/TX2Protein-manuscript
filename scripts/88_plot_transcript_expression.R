
# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# Load data ---------------------------------------------------------------

expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))

# Main --------------------------------------------------------------------

# Only include genes expressed across all samples
genes_to_include <- 
  expression_gene_set %>% 
  dplyr::select(annot_gene_id, sample) %>% 
  group_by(annot_gene_id) %>%
  summarize(n_samples = n_distinct(sample)) %>%
  dplyr::filter(n_samples == 9)

# Only include transcripts expressed across 2 or more samples
transcripts_to_include <-
  expression_gene_set %>% 
  dplyr::filter(annot_gene_id %in% genes_to_include$annot_gene_id,
                transcript_novelty %in% c("NNC", "NIC"))

# Plot number of samples each transcript is expressed in
no_samples_transcript_expressed_plot <-
  expression_gene_set %>% 
  dplyr::filter(unique_id %in% transcripts_to_include$unique_id) %>% 
  dplyr::select(unique_id, sample) %>%
  group_by(unique_id) %>%
  summarize(n_samples = n_distinct(sample)) %>%
  group_by(n_samples) %>% 
  summarize(n_transcripts = dplyr::n()) %>% 
  ggplot(aes(x = n_samples,
             y = n_transcripts)) + 
  geom_col(colour = "black",
           fill = "lightgrey") +
  geom_text(aes(label = n_transcripts), 
            position = position_dodge(width=0.9), vjust= -0.6, size = 4) +
  scale_x_continuous(breaks = c(1:9)) +
  labs(x = "",
       y = "Number of genes") +
  theme_light() +
  theme(axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# Plot mean expression of each transcript

# get number of samples expressing each transcript
sample_expression <-
  expression_gene_set %>% 
  dplyr::filter(unique_id %in% transcripts_to_include$unique_id) %>% 
  dplyr::select(unique_id, sample) %>%
  group_by(unique_id) %>%
  summarize(n_samples = n_distinct(sample))

# plot
expression_per_transcript_plot <-
  expression_gene_set %>% 
  dplyr::filter(unique_id %in% transcripts_to_include$unique_id) %>% 
  dplyr::select(unique_id, sample, RPM) %>% 
  aggregate(data = ., RPM ~ unique_id, FUN = mean) %>% 
  dplyr::left_join(., sample_expression, by = join_by(unique_id)) %>% 
  ggplot(aes(x = n_samples,
             y = RPM, 
             group = n_samples)) + 
  geom_boxplot(colour = "black", fill = "lightgrey", show.legend = F) +
  scale_x_continuous(breaks = c(1:9)) +
  labs(x = "No. samples expressed in (total = 9)",
       y = "Expression (RPM)") +
  ylim(0, 20) +
  theme_light() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# Arrange plots
transcript_expression_plot <-
  ggpubr::ggarrange(no_samples_transcript_expressed_plot,
                    expression_per_transcript_plot,
                    ncol = 1, 
                    align = "v",
                    labels = c("A", "B"), 
                    font.label = list(size = 20))

# Save data ---------------------------------------------------------------

ggsave(
  plot = transcript_expression_plot, 
  filename = "SuppFig2_transcript_expression_plot.png", 
  path = here::here("results", "plots"), 
  width = 8, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)
