# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# Load data ---------------------------------------------------------------

mass_Spec_hits <- read.table(here::here("results", "mass_spec_hits.txt"), col.names = c("Transcript", "no_use", "number_of_occurrences"), sep = "\t")

expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))

druggable_genome <- read.csv(here::here("raw_data", "pharos_data_download", "query results.csv"))

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


# Number of transcripts per transcript category
transcripts_to_include <- 
  expression_gene_set %>% 
  dplyr::select(annot_gene_id, unique_id, transcript_novelty, sample, reads) %>% 
  dplyr::filter(annot_gene_id %in% genes_to_include$annot_gene_id,
                transcript_novelty %in% c("NNC", "NIC"),
                reads > 1) %>% 
  dplyr::select(unique_id, sample) %>% 
  group_by(unique_id) %>%
  summarize(n_samples = n_distinct(sample)) %>% 
  dplyr::filter(n_samples > 1)




# filter data
final_transcript_set_MS <-
  mass_Spec_hits %>% 
  dplyr::filter(Transcript %in% transcripts_to_include$unique_id)

plot_MS_support <- 
  final_transcript_set_MS %>% 
  group_by(number_of_occurrences) %>%
  summarize(n_transcripts = n_distinct(Transcript)) %>% 
  
  ggplot(aes(x = number_of_occurrences,
             y = n_transcripts, 
             group = number_of_occurrences)) + 
  geom_col(aes(fill = factor(ifelse(number_of_occurrences > 0, "green", "lightgrey"))),
           colour = "black") +
  geom_text(aes(label = n_transcripts), 
            position = position_dodge(width=0.9), vjust= -0.6, size = 4) +
  labs(x = "No. mass spec data sets (total = 30)",
       y = "No. transcripts") +
  scale_fill_manual(values = c("green" = "#00B358", "lightgrey" = "lightgrey")) +
  theme_light() +
  guides(fill = "none") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

#### plot genes with MS support
transcript_set_MS_no_0 <- final_transcript_set_MS %>% dplyr::filter(number_of_occurrences > 0) 

gene_transcript <- 
  expression_gene_set %>% dplyr::select(annot_gene_name, RPM, unique_id) %>% 
  dplyr::filter(unique_id %in% transcripts_to_include$unique_id) %>% 
  aggregate(data = ., RPM ~ unique_id + annot_gene_name, FUN = mean)

final_MS_gene <- transcript_set_MS_no_0 %>% dplyr::left_join(., gene_transcript, by = c("Transcript" = "unique_id")) %>% 
  aggregate(data = ., RPM ~ annot_gene_name, FUN = sum) %>% 
  arrange(RPM)

plot_gene_novelty <- 
  final_MS_gene %>% 
  mutate(category = ifelse(RPM > 20, "Above 20 RPM", "Below or Equal 20 RPM")) %>% 
  ggplot(
    aes(x = factor(annot_gene_name, levels = annot_gene_name),
        y = RPM)) +
  geom_col(colour = "black",
           fill = "lightgrey") +
  labs(x = "Gene",
       y = "Expression of novel protein coding transcripts (RPM)") +
  coord_flip() +
  facet_wrap(~category, scales = "free", nrow = 1) +
  theme_light() +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0.0)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7))

########### How much does this novel protein coding transcript correspond to in relative terms of expression ##########
genes_with_novel_prot <- 
  expression_gene_set %>% 
  dplyr::filter(unique_id %in% transcript_set_MS_no_0$Transcript) %>% 
  dplyr::select(annot_gene_name) %>% 
  unique()

total_expression_per_sample <-
  expression_gene_set %>% 
  dplyr::filter(annot_gene_name %in% genes_with_novel_prot$annot_gene_name) %>% 
  aggregate(data = ., RPM ~ annot_gene_name + sample, FUN = sum)


novel_expression <-
  expression_gene_set %>% 
  dplyr::filter(unique_id %in% transcript_set_MS_no_0$Transcript) %>% 
  aggregate(data = ., RPM ~ annot_gene_name + sample, FUN = sum)

combined <- 
  total_expression_per_sample %>% 
  dplyr::left_join(., novel_expression, by = join_by(annot_gene_name, sample)) %>% 
  replace(., is.na(.), 0) %>% 
  dplyr::mutate(novel_expression_perc = (RPM.y / RPM.x) *100)

top_genes <- 
  combined %>% 
  aggregate(data = ., novel_expression_perc ~ annot_gene_name, FUN = mean) %>% 
  arrange(desc(novel_expression_perc)) %>% 
  slice_max(novel_expression_perc, n = 20)

plot_novel_expression_all <-
  combined %>% 
  ggplot(
    aes(x = reorder(annot_gene_name, -novel_expression_perc),
        y = novel_expression_perc,
        group = annot_gene_name)) +
  geom_boxplot(colour = "black", fill = "lightgrey", show.legend = F, outlier.shape = NA) +
  labs(x = "Gene",
       y = "Relative expression of novel protein coding transcripts (%)") +
  theme_light() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1)
)

plot_novel_expression_top <-
  combined %>% 
  dplyr::filter(annot_gene_name %in% top_genes$annot_gene_name) %>% 
  ggplot(
    aes(x = reorder(annot_gene_name, -novel_expression_perc),
        y = novel_expression_perc,
        group = annot_gene_name)) +
  geom_boxplot(colour = "black", fill = "lightgrey", show.legend = F, outlier.shape = NA) +
  geom_point() +
  labs(x = "Gene",
       y = "Relative expression of novel protein coding transcripts (%)") +
  theme_light() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1)
  )

#####

plot_novel_expression_with_zoom <-
  combined %>% 
  ggplot(
    aes(x = as.numeric(reorder(annot_gene_name, -novel_expression_perc)),
        y = novel_expression_perc,
        group = annot_gene_name)) +
  geom_boxplot(colour = "black", fill = "lightgrey", show.legend = F, outlier.shape = NA) +
  scale_x_continuous(
    breaks = 1:length(levels(combined$annot_gene_name)),
    label = levels(combined$annot_gene_name)) +
  labs(x = "Gene",
       y = "Relative expression of novel protein coding transcripts (%)") +
  facet_zoom(x = annot_gene_name %in% top_genes$annot_gene_name) +
  theme_light() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1)
  )

#####


# Save data ---------------------------------------------------------------

ggsave(
  plot = plot_MS_support, 
  filename = "plot_MS_support.png", 
  path = here::here("results", "plots"), 
  width = 8, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)

ggsave(
  plot = plot_gene_novelty, 
  filename = "plot_gene_novelty.png", 
  path = here::here("results", "plots"), 
  width = 16, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)

ggsave(
  plot = plot_novel_expression_all, 
  filename = "plot_novel_expression_all.png", 
  path = here::here("results", "plots"), 
  width = 16, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)

ggsave(
  plot = plot_novel_expression_top, 
  filename = "plot_novel_expression_top.png", 
  path = here::here("results", "plots"), 
  width = 12, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)