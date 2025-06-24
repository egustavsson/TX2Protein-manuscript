# Load libraries ----------------------------------------------------------
library(ORFik)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(ggtranscript)

# Arguments ---------------------------------------------------------------
gene_of_interest <- "ALG3"
MANE <- "ENST00000397676"
fa_path <- "/home/MinaRyten/Emil/references/GRCh38.primary_assembly.genome.fa"
rds_path <- here::here("results", "unique_transcripts.rds")
TX2P_path <- here::here("results", "TX2P_output", "EPI_Revised_main.xlsx")
expression_path <- here::here("results", "expression_gene_set.rds")

# Load data ---------------------------------------------------------------

# full transcript model with gene_name metadata
all_features <- readRDS(rds_path)
exon_gr <- all_features %>% dplyr::filter(type == "exon") %>% GRanges()
TX2P_count <- read_excel(TX2P_path, sheet = "Count_Found_in_mass_Datasets", col_names = c("Transcript_ID", "occurences"))
expression <- readRDS(here::here("results", "expression_gene_set.rds"))

# also download and load reference annotation 
ref_path <- here::here(tempdir(), "gencode.v38.annotation.gtf.gz")

if(!file.exists(ref_path)) {
  
  download.file(
    url = paste0(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/", 
      "gencode.v38.annotation.gtf.gz"
    ),
    destfile = ref_path
  )
  
}

R.utils::gunzip(ref_path, remove = TRUE)

ref <- rtracklayer::import(stringr::str_remove(ref_path, "\\.gz"))


# Functions ---------------------------------------------------------------

get_lr_tx_of_interest <- function(lr, pb_ids) {
  
  lr <- lr[lr$transcript_id %in% pb_ids]
  lr_exons <- lr[lr$type == "exon"]
  lr_cds <- lr[lr$type == "CDS"]
  
  # return as list as we need both exons and cds
  lr_exons_cds <- list(
    exons = lr_exons, 
    cds = lr_cds
  )
  
  return(lr_exons_cds)
  
}



get_mane <- function(ref, mane_id) {
  
  # remove any NA transcript ids (i.e. type == "gene")
  mane <- ref[!is.na(ref$transcript_id)] 
  
  # remove to .XX after ENST
  GenomicRanges::mcols(mane)[["transcript_id"]] <- 
    GenomicRanges::mcols(mane)[["transcript_id"]] %>% 
    stringr::str_remove("\\..*")
  
  mane <- mane[mane$transcript_id == mane_id, ]
  mane_exons <- mane[mane$type == "exon"]
  mane_cds <- mane[mane$type == "CDS"]
  
  mane_exons_cds <- list(
    exons = mane_exons, 
    cds = mane_cds
  )
  
  return(mane_exons_cds)
  
}



# annotate CDS
annotate_longest_orf_cds <- function(exon_gr,
                                     gene_name,
                                     fasta_path,
                                     tx_id_col = "transcript_id",
                                     gene_name_col = "gene_name",
                                     keep_only_with_orf = TRUE) {
  # Subset exons for gene of interest
  gene_exons <- exon_gr[mcols(exon_gr)[[gene_name_col]] == gene_name]
  if (length(gene_exons) == 0) return(gene_exons)
  
  # Build GRangesList per transcript
  grl <- split(gene_exons, mcols(gene_exons)[[tx_id_col]])
  grl <- ORFik:::sortPerGroup(grl)
  
  # Extract transcript sequences
  seqs <- ORFik:::txSeqsFromFa(grl, Rsamtools::FaFile(fasta_path))
  
  # Predict ORFs, allow multiple per transcript
  cds_grl <- findMapORFs(
    grl = grl,
    seqs = seqs,
    startCodon = "ATG",
    longestORF = TRUE,
    groupByTx = TRUE
  )
  
  # Step 1: Compute ORF widths per `mcols(x)$names` (e.g. transcript_1, transcript_2)
  orf_widths_by_name <- lapply(cds_grl, function(x) {
    name_vec <- mcols(x)$names
    unique_names <- unique(name_vec)
    
    sapply(unique_names, function(nm) {
      sum(width(x[name_vec == nm]))
    })
  })
  
  # Step 2: Flatten to data.frame and keep longest ORF per transcript
  orf_df <- bind_rows(
    lapply(orf_widths_by_name, function(x) {
      longest <- which.max(x)
      tibble(
        orf_name = names(x)[longest],
        orf_width = x[longest]
      )
    }),
    .id = "transcript_index"
  )
  
  # Step 3: Flatten and annotate ORFs
  cds_flat <- unlist(cds_grl, use.names = FALSE)
  mcols(cds_flat)$type <- "CDS"
  mcols(cds_flat)$gene_name <- gene_name
  mcols(cds_flat)$transcript_id <- sub("_\\d+$", "", mcols(cds_flat)$names)
  
  # Step 4: Filter to longest ORF per transcript
  cds_flat_final <- cds_flat[mcols(cds_flat)$names %in% orf_df$orf_name, ]
  
  # Optionally filter exons to only those with ORFs
  if (keep_only_with_orf) {
    gene_exons <- gene_exons[mcols(gene_exons)[[tx_id_col]] %in% unique(mcols(cds_flat_final)$transcript_id)]
  }
  
  return(c(gene_exons, cds_flat_final))
}


# plot transcripts
plot_diff <- function(lr_exons_cds, 
                      mane_exons_cds, 
                      lr_mane_diffs, 
                      MANE_canonical = "MANE"
) {
  
  # merge mane and lr data and convert to data.frame() for plotting
  # convert transcript_id to factor to make sure mane is at top
  transcript_order <- c(
    lr_exons_cds$exons$transcript_id %>% unique(),
    mane_exons_cds$exons$transcript_id %>% unique()
  )
  
  lr_mane_exons_df <- c(lr_exons_cds$exons, mane_exons_cds$exons) %>% 
    data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  lr_mane_cds_df <- c(lr_exons_cds$cds, mane_exons_cds$cds) %>% 
    data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  
  # plot diff plot
  diff_plot <- lr_mane_exons_df %>% 
    ggplot(aes(
      xstart = start, 
      xend = end, 
      y = transcript_id
    )) + 
    geom_range(
      height = 0.25, 
      fill = "white"
    ) + 
    geom_range(
      data = lr_mane_cds_df
    ) + 
    geom_intron(
      data = to_intron(lr_mane_exons_df, "transcript_id"), 
      aes(strand = strand), 
      arrow.min.intron.length = 400
    ) + 
    geom_range(
      data = lr_mane_diffs, 
      aes(
        fill = diff_type, 
        colour = diff_type
      ), 
      alpha = 0.2, 
      linetype = 2
    ) + 
    scale_y_discrete(name = "Transcript ID") + 
    scale_x_continuous(name = "Genomic position") + 
    scale_fill_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    scale_colour_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    theme_bw() + 
    theme(legend.position = "top")
  
  return(diff_plot)
  
}


# Main --------------------------------------------------------------------

mane_exons_cds <- get_mane(ref, mane_id = MANE)

# only include transcripts with MS support
exon_gr_final <- exon_gr[exon_gr$transcript_id %in% TX2P_count$Transcript_ID,]


result <- annotate_longest_orf_cds(exon_gr = exon_gr_final,
                                   gene_name = gene_of_interest,
                                   fasta_path = fa_path)


lr_exons_cds <- get_lr_tx_of_interest(
  lr = result, 
  pb_ids = result$transcript_id
  )


# obtain differences between MANE and lr exons
lr_mane_diffs <- 
  ggtranscript::to_diff(
    exons = lr_exons_cds$exons %>% data.frame(), 
    ref_exons = mane_exons_cds$exons %>% data.frame(), 
    group_var = "transcript_id"
  )

lr_mane_diff_plot <- 
  plot_diff(
    lr_exons_cds = lr_exons_cds, 
    mane_exons_cds = mane_exons_cds, 
    lr_mane_diffs = lr_mane_diffs
  )

# process expression
expression_relative <- expression %>%
  dplyr::select(!c(gene_ID,
                   annot_gene_id, 
                   n_exons, 
                   length, 
                   gene_novelty, 
                   transcript_novelty, 
                   ISM_subtype, 
                   status, 
                   tissue, 
                   age, 
                   sex)) %>% 
  group_by(sample, annot_gene_name) %>%
  mutate(RPM_relative = 100 * RPM / sum(RPM)) %>%
  ungroup()

# plot expression
expression_relative %>%
  mutate(transcript_id = stringr::str_remove(unique_id, "\\..*")) %>%
  filter(transcript_id %in% c(lr_mane_diffs$transcript_id, MANE)) %>% 
  ggplot(aes(unique_id, y = RPM_relative)) +
  geom_boxplot(
    width = 0.5, outlier.shape = NA
  ) +
  geom_point() +
  coord_flip() +
  labs(x = "", 
       y = "Relative expression per transcript") +
  scale_y_continuous(labels = function(x) paste0(x, '%')) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

# Create shared transcript_id factor
transcript_order <- c(
  lr_exons_cds$exons$transcript_id %>% unique(),
  mane_exons_cds$exons$transcript_id %>% unique()
)

# Prep expression data for boxplot
expression_for_plot <- expression_relative %>%
  mutate(transcript_id = stringr::str_remove(unique_id, "\\..*")) %>%
  filter(transcript_id %in% transcript_order) %>%
  mutate(transcript_id = factor(transcript_id, levels = transcript_order))

# Expression boxplot (rotated to match ggtranscript plot)
expr_boxplot <- ggplot(expression_for_plot, aes(y = transcript_id, x = RPM_relative)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0, height = 0.2, alpha = 0.5) +
  scale_x_continuous(name = "Relative expression (%)", labels = scales::percent_format(scale = 1)) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),  # hide y-labels to match ggtranscript
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
  )

# Combine plots
combined_plot <- lr_mane_diff_plot + expr_boxplot + plot_layout(widths = c(3, 1))

ggsave(
  plot = combined_plot, 
  filename = "04a_transcript_expression_plot.svg", 
  path = here::here("results", "plots"), 
  width = 12, 
  height = 5, 
  dpi = 600, 
  bg = "white"
)
