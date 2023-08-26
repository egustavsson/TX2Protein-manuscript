# Make a histogram as a stacked bar chart with different transcript categories


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

transcript_category_plot <-
  expression_gene_set %>% 
  dplyr::select(annot_gene_name, transcript_novelty, unique_id) %>% 
  distinct() %>% 
  group_by(annot_gene_name, transcript_novelty) %>%
  summarize(count = n_distinct(unique_id)) %>% 
  dplyr::filter(transcript_novelty != "Genomic") %>% 

  # Plot histogram
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
        legend.background = element_rect(fill = alpha("white", 0.0)),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))


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
