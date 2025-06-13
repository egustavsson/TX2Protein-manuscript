# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggpubr)
})

# Load data ---------------------------------------------------------------

expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))

# Main --------------------------------------------------------------------

# Plot the number of samples in which each gene is expressed  
no_samples_expressed_plot <-
  expression_gene_set %>% 
  dplyr::select(annot_gene_name, sample) %>% 
  group_by(annot_gene_name) %>%
  summarize(n_samples = n_distinct(sample)) %>%
  group_by(n_samples) %>% 
  summarize(n_genes = dplyr::n()) %>% 
  ggplot(aes(x = n_samples,
             y = n_genes)) + 
  geom_col(colour = "black",
           fill = "lightgrey") +
  geom_text(aes(label = n_genes), 
            position = position_dodge(width=0.9), vjust= -0.6, size = 6) +
  scale_x_continuous(breaks = c(1:9)) +
  geom_vline(xintercept = 8.5, linetype = "dashed", color = "black", size = 0.7) + # Add dashed line
  labs(x = "",
       y = "Number of genes") +
  ylim(0, 520) +
  theme_light() +
  theme(axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# get number of samples each gene is expressed in
samples_per_gene <-
  expression_gene_set %>% 
  dplyr::select(annot_gene_name, annot_gene_id, sample) %>% 
  group_by(annot_gene_name, annot_gene_id) %>%
  summarize(n_samples = n_distinct(sample))

# Plot mean expression of each gene
expression_per_gene_plot <-
  expression_gene_set %>% 
  dplyr::select(annot_gene_id, sample, RPM) %>% 
  aggregate(data = ., RPM ~ annot_gene_id + sample, FUN = sum) %>% 
  dplyr::left_join(., samples_per_gene, by = "annot_gene_id") %>% 
  ggplot(aes(x = n_samples,
             y = RPM, 
             group = n_samples)) + 
  geom_boxplot(colour = "black", fill = "lightgrey", show.legend = F, outlier.shape = NA) +
  scale_x_continuous(breaks = c(1:9)) +
  geom_vline(xintercept = 8.5, linetype = "dashed", color = "black", size = 0.7) + # Add dashed line
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
  plot = gene_expression_plot, 
  filename = "SFig1_gene_expression_plot.png", 
  path = here::here("results", "plots"), 
  width = 8, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)
