# This script is used to filter the TX2P output
# 
# We ran all the transcripts through TX2P but will only include
# those expressed across all samples included and with a read cut-off

# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
  library(ggrepel)
  library(ggforce)
  library(ggrepel)
})

# Arguments ---------------------------------------------------------------

args <- list(
  path_to_TX2P = here::here("results", "TX2P_output", "EPI_Revised_main.xlsx"),
  path_to_expression = here::here("results", "expression_gene_set.rds")
)

# Load data ---------------------------------------------------------------

# TX2P output data
# This might need updating depending on the final output format we decide on for TX2P
# TX2P_count contains the number of times each transcript had a positive peptide detection
# TX2P_merged has all the data
TX2P_out <- list(
  TX2P_count = read_excel(args$path_to_TX2P, sheet = "Count_Found_in_mass_Datasets", col_names = c("Transcript_ID", "occurences"))
  )

# This data is generated through 01b_process_ENCODE_expression.R
expression_gene_set <- readRDS(here::here("results", "expression_gene_set.rds"))


# Functions ---------------------------------------------------------------

TXfilter <- function(data, gene_expression_sample_count, min_read_count, allowed_novelty_types, min_transcript_occurrences) {
  
  # Step 1: Identify genes expressed in at least gene_expression_sample_count samples
  selected_genes <- 
    data %>%
    dplyr::select(annot_gene_id, sample) %>%
    group_by(annot_gene_id) %>%
    summarize(sample_count = n_distinct(sample), .groups = "drop") %>%
    dplyr::filter(sample_count >= gene_expression_sample_count) %>%
    dplyr::pull(annot_gene_id)
  
  # Step 2: Identify transcripts found in at least min_transcript_occurrences samples
  selected_transcripts <- 
    data %>%
    dplyr::select(unique_id, sample) %>%
    group_by(unique_id) %>%
    summarize(transcript_occurrences = n_distinct(sample), .groups = "drop") %>%
    dplyr::filter(transcript_occurrences >= min_transcript_occurrences) %>%
    dplyr::pull(unique_id)
  
  # Step 3: Filter transcripts of interest
  filtered_transcripts <- 
    data %>%
    dplyr::filter(
      annot_gene_id %in% selected_genes,
      transcript_novelty %in% allowed_novelty_types, 
      reads >= min_read_count,
      unique_id %in% selected_transcripts
    )
  
  return(filtered_transcripts)
}

# Main --------------------------------------------------------------------

filtered_data <- TXfilter(
  data = expression_gene_set,
  gene_expression_sample_count = 9, 
  min_transcript_occurrences = 1,
  min_read_count = 2,
  allowed_novelty_types = c("NNC", "NIC")
)


# get filtered data and annotate with MS support
filtered_data_w_MS_support <- filtered_data %>% 
  dplyr::left_join(
    TX2P_out$TX2P_count %>% 
      dplyr::select(Transcript_ID, occurences) %>% 
      distinct(), 
    by = c("unique_id" = "Transcript_ID")
  ) %>% 
  dplyr::mutate(occurences = tidyr::replace_na(occurences, 0))

# These plots are to show how many transcripts that are supported by MS data and
# by how many MS datasets they are supported by

# --------------------------------------------
# plot support by number of MS datasets
# --------------------------------------------
transcripts_number_of_MS_data_support_plot <-
  filtered_data_w_MS_support %>% 
  dplyr::select(unique_id, occurences) %>% 
  unique() %>% 
  group_by(occurences) %>%
  summarize(n_transcripts = n_distinct(unique_id)) %>% 
  
  ggplot(aes(x = occurences,
             y = n_transcripts, 
             group = occurences)) + 
  geom_col(aes(fill = factor(ifelse(occurences > 0, "blue", "lightgrey"))),
           colour = "black") +
  geom_text(aes(label = n_transcripts), 
            position = position_dodge(width=0.9), vjust= -0.6, size = 4) +
  labs(x = "No. mass spec data sets (total = 27)",
       y = "No. transcripts") +
  scale_fill_manual(values = c("blue" = "#8ed6ffff", "lightgrey" = "#d7d7d780")) +
  theme_light() +
  guides(fill = "none") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# --------------------------------------------
# plot support by MS as a binary value
# --------------------------------------------
transcripts_MS_support_plot <-
  filtered_data_w_MS_support %>% 
  dplyr::select(unique_id, occurences) %>% 
  dplyr::mutate(MS_support = occurences > 0) %>%  # Ensure logical type
  distinct() %>% 
  group_by(MS_support) %>%
  summarize(n_transcripts = n_distinct(unique_id), .groups = "drop") %>%
  mutate(
    perc = 100 * n_transcripts / sum(n_transcripts),                
    label = paste0(n_transcripts, " (", sprintf("%.1f", perc), "%)"),  
    fill_color = ifelse(MS_support, "#8ed6ffff", "#d7d7d780")
  ) %>%
  
  ggplot(aes(x = MS_support, y = n_transcripts, fill = fill_color)) + 
  geom_col(colour = "black", width = 0.6) + 
  geom_text(aes(label = label), vjust = -0.6, size = 5) +
  labs(title = "Transcripts with MS support",
       x = "Supported by mass spectrometry data", 
       y = "No. transcripts") +
  scale_fill_identity() +  # Use identity scale to apply pre-defined colors
  theme_bw() +
  guides(fill = "none") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5)
  )


# --------------------------------------------
# Plot cummulative MS support
# --------------------------------------------

# Aggregate RPM values
aggregated_data <- filtered_data_w_MS_support %>%
  aggregate(data = ., RPM ~ unique_id + occurences, FUN = mean) %>% 
  # Create a flag for occurrences > 0
  mutate(occurence_flag = ifelse(occurences > 0, 1, 0))

# Define RPM thresholds
threshs <- seq(from = 0, to = ceiling(max(aggregated_data$RPM)), by = 0.1)

# Create a data frame to store results
res <- data.frame()

# Get total number of transcripts (denominator) and max proportion
total_transcripts <- nrow(aggregated_data)
max_proportion <- sum(aggregated_data$occurence_flag) / total_transcripts  # Proportion of transcripts with occurence_flag == 1

# Loop through each threshold and calculate the cumulative proportion
for(i in threshs) {
  threshold_data <- aggregated_data %>%
    dplyr::filter(RPM <= i) %>%  # Cumulative count of transcripts with RPM ≤ threshold
    summarise(RPM_cutoff = i, cumulative_proportion = sum(occurence_flag) / total_transcripts)  # Normalize by total transcripts
  
  # Combine results for each threshold
  res <- rbind(res, threshold_data)
}

# Plotting the results
cumulative_MS_support_plot <-
  res %>% 
  ggplot(aes(x = RPM_cutoff, y = cumulative_proportion)) +
  geom_area(fill = "#B0E2FF",
            alpha = 0.5,
            color = "#607B8B",
            lwd = 1,
            linetype = 1) +
  ggtitle("Cumulative Proportion of Transcripts with Mass Spec Support by RPM") +
  scale_x_continuous(
    name = "Mean expression (RPM)"
  ) +
  scale_y_continuous(
    name = "Cumulative Proportion of Transcripts with Mass Spec Support",
    limits = c(0, max_proportion),  
    breaks = seq(0, max_proportion, by = 0.05),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_zoom(xlim = c(0, 5)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12)
  )


# --------------------------------------------
# Plot expression by gene
# --------------------------------------------

# These plots shows examples of genes with a large proportion of transcription being
# novel and supported by mass spec

# generate total gene expression for each gene.
# For denominator i use all transcripts regardless of MS support and for numerator i use novel with MS support

# total expression
total_gene_expression <- expression_gene_set %>%
  dplyr::select(annot_gene_id, annot_gene_name, sample, RPM) %>% 
  aggregate(data =., RPM ~ annot_gene_id + annot_gene_name + sample, FUN = sum)

# novel tx expression, only including those with MS support
novel_tx_expression <- filtered_data_w_MS_support %>% 
  dplyr::filter(occurences > 0) %>%
  dplyr::select(annot_gene_id, annot_gene_name, sample, RPM) %>% 
  aggregate(data =., RPM ~ annot_gene_id + annot_gene_name + sample, FUN = sum)

# calculate proportion novel expression per sample
prop_novel_expression <- 
  total_gene_expression %>% 
  dplyr::filter(annot_gene_id %in% novel_tx_expression$annot_gene_id) %>% # this is that i only calculate for genes with novel expression
  dplyr::left_join(., novel_tx_expression, by = c("annot_gene_name" = "annot_gene_name", "sample" = "sample")) %>% 
  replace(., is.na(.), 0) %>% 
  dplyr::mutate(novel_expression_perc = (RPM.y / RPM.x) *100)

# create data frame of percentage novel expression
summary_expression_df <- prop_novel_expression %>%
  group_by(annot_gene_id.x, annot_gene_name) %>%
  summarise(
    mean_novel_expression_perc = mean(novel_expression_perc, na.rm = TRUE),
    sd_novel_expression_perc = sd(novel_expression_perc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(
    annot_gene_id = annot_gene_id.x)


# highlight the top 5 genes
highlight_genes <- prop_novel_expression %>%
  arrange(desc(novel_expression_perc)) %>%
  distinct(annot_gene_name, .keep_all = TRUE) %>% 
  slice_head(n = 5) %>%
  pull(annot_gene_name)
  
  
# plot all genes
prop_novel_expression_all_plot <-
  prop_novel_expression %>% 
  aggregate(data = ., novel_expression_perc ~ annot_gene_name, FUN = mean) %>% 
  arrange(desc(novel_expression_perc)) %>%
  ggplot(
    aes(x = reorder(annot_gene_name, -novel_expression_perc),
        y = novel_expression_perc)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.6) +
  geom_point(size = 2, shape = 21, color = "black", fill = "#8ed6ffff") + # Darker fill
  geom_label_repel(aes(label = ifelse(annot_gene_name %in% highlight_genes, annot_gene_name, NA)),
                   size = 4, box.padding = 0.3, point.padding = 0.2, 
                   label.size = 0.3, fill = "#8ed6ffff", max.overlaps = Inf) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_bw() +
  labs(x = "", 
       y = "Mean Relative Expression of \n Novel Protein Coding Transcript") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# plot top 5 genes
prop_novel_expression_top_plot <-
  prop_novel_expression %>% 
  dplyr::filter(annot_gene_name %in% highlight_genes) %>%
  ggplot(
    aes(x = reorder(annot_gene_name, -novel_expression_perc),
        y = novel_expression_perc,
        group = annot_gene_name)) +  
  geom_boxplot(show.legend = FALSE, 
               outlier.shape = NA,
               na.rm = TRUE,
               width = 0.6,
               fill = "#8ed6ffff") +  
  geom_point(show.legend = FALSE,
             na.rm = TRUE,
             size = 2,
             alpha = 0.6,
             shape = 21,
             fill = "#8ed6ffff") +  
  labs(
    x = "Gene",
    y = "Mean Relative Expression of \n Novel Protein Coding Transcript") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold.italic", hjust = 1, angle = 45),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5)
  )

# Arrange plots with a single title centered
combined_plot <- ggpubr::ggarrange(
  plotlist = list(prop_novel_expression_all_plot, prop_novel_expression_top_plot), 
  nrow = 1, 
  align = "h", 
  widths = c(1.5, 1)
)

# Add the title to the combined plot using annotate_figure
prop_novel_expression_plot <- ggpubr::annotate_figure(
  combined_plot, 
  top = text_grob("Novel Protein Coding Transcript Expression by Gene", 
                  face = "bold", size = 20, hjust = 0.5)
)

# Save data ---------------------------------------------------------------

# plot support by MS as a binary value
ggsave(
  plot = transcripts_MS_support_plot, 
  filename = "04a_transcripts_MS_support_plot.svg", 
  path = here::here("results", "plots"), 
  width = 6, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)

# plot support by number of MS data sets (supplementary)
ggsave(
  plot = transcripts_number_of_MS_data_support_plot, 
  filename = "04a_transcripts_number_of_MS_data_support_plot.svg", 
  path = here::here("results", "plots"), 
  width = 10, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)

# plot MS support cummulatively
ggsave(
  plot = cumulative_MS_support_plot, 
  filename = "04a_cumulative_MS_support_plot.svg", 
  path = here::here("results", "plots"), 
  width = 11, 
  height = 9, 
  dpi = 600, 
  bg = "white"
)

# plot novel protein coding by gene
ggsave(
  plot = prop_novel_expression_plot, 
  filename = "04a_prop_novel_expression_plot.svg", 
  path = here::here("results", "plots"), 
  width = 18, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)

# table with novel expression per sample
summary_expression_df %>%
  write.csv(here("results", "novel_expression.csv"), row.names = FALSE)




