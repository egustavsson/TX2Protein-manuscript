
# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# Load data ---------------------------------------------------------------

expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))

# Main --------------------------------------------------------------------

# fill colour to use
fill_colour <- c("Known" = "#4d9221",
                 "ISM" = "#74add1",
                 "NIC" = "#d53e4f",
                 "NNC" = "#b2abd2")

# Only include genes expressed across all samples
genes_to_include <- 
  expression_gene_set %>% 
  dplyr::select(annot_gene_name, annot_gene_id, sample, RPM) %>% 
  group_by(annot_gene_id) %>%
  summarize(n_samples = n_distinct(sample)) %>%
  dplyr::filter(n_samples == 9)
  
# Number of transcripts per transcript category
transcript_category_data <-
  expression_gene_set %>% 
  dplyr::select(annot_gene_id, annot_gene_name, transcript_novelty, unique_id, reads) %>% 
  distinct() %>% 
  dplyr::filter(transcript_novelty != "Genomic",
                annot_gene_id %in% genes_to_include$annot_gene_id,
                reads > 1) %>% 
  group_by(annot_gene_id, transcript_novelty) %>%
  summarize(count = n_distinct(unique_id))


# Plot total transcript distribution
transcript_total_plot <-
  transcript_category_data %>% 
  aggregate(data = ., count ~ transcript_novelty, FUN = sum) %>% 
  ggplot(aes(
    y = count, 
    x = factor(transcript_novelty, 
                  levels = c("Known", "ISM", "NIC", "NNC")),
    fill = transcript_novelty)) + 
  geom_col(colour = "black") +
  geom_text(aes(label = count), 
            position = position_dodge(width=0.9), vjust= 2, size = 8) +
  labs(x = "Genes",
       y = "Transcripts per category",
       fill = "Transcript category") +
  scale_fill_manual(values = fill_colour) +
  theme_light() +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0.0)),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))


# Plot per gene in a histogram
transcript_category_per_gene_plot <-
  transcript_category_data %>% 
  ggplot(aes(
    x = count, 
    fill = factor(transcript_novelty, 
                  levels = c("Known", "ISM", "NIC", "NNC")))) + 
  geom_histogram(colour = "black",
                 bins = 100) +
  labs(x = "Genes",
       y = "Distribution of transcripts per category",
       fill = "Transcript category") +
  scale_fill_manual(values = fill_colour) +
  theme_light() +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0.0)),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))


#### This is to plot number of samples each transcript is expressed in ####

transcripts_to_include <-
  expression_gene_set %>% 
  dplyr::filter(transcript_novelty %in% c("NNC", "NIC"),
                annot_gene_id %in% genes_to_include$annot_gene_id) %>% 
  dplyr::select(unique_id, sample) %>% 
  group_by(unique_id) %>%
  summarize(n_samples = n_distinct(sample)) %>%
  dplyr::filter(n_samples > 1)


# Plot how many samples is each transcript is expressed in
no_samples_transcript_expressed_plot <-
  expression_gene_set %>% 
  #dplyr::filter(transcript_novelty %in% c("NNC", "NIC"),
 #               annot_gene_id %in% genes_to_include$annot_gene_id) %>% 
  dplyr::filter() %>% 
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
       y = "Number of transcripts") +
  #ylim(0, 520) +
  scale_fill_manual(values = fill_colour) +
  theme_light() +
  theme(axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# Plot mean expression of each gene
expression_per_gene_plot <-
  expression_gene_set %>% 
  dplyr::select(unique_id, sample, RPM) %>% 
  aggregate(data = ., RPM ~ unique_id + sample, FUN = sum) %>% 
  dplyr::left_join(., sample_expression, by = join_by(annot_gene_id)) %>% 
  ggplot(aes(x = n_samples,
             y = RPM, 
             group = n_samples)) + 
  geom_boxplot(colour = "black", fill = "lightgrey", show.legend = F, outlier.shape = NA) +
  ggpubr::stat_compare_means(paired = F) +
  scale_x_continuous(breaks = c(1:9)) +
  labs(x = "No. samples expressed in (total = 9)",
       y = "Expression (RPM)") +
  ylim(0, 100) +
  theme_light() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# Arrange plots
gene_expression_plot <-
  ggpubr::ggarrange(no_samples_expressed_plot,
                    expression_per_gene_plot,
                    ncol = 1, 
                    align = "v",
                    labels = c("A", "B"), 
                    font.label = list(size = 20))

# Save data ---------------------------------------------------------------

ggsave(
  plot = transcript_category_plot, 
  filename = "transcript_category_plot.png", 
  path = here::here("results", "plots"), 
  width = 8, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)

ggsave(
  plot = transcript_category_per_gene_plot, 
  filename = "transcript_category_per_gene_plot.png", 
  path = here::here("results", "plots"), 
  width = 8, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)
