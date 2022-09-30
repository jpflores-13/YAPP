## ATAC peaks at static, gained, and lost anchors with DESeq default normalization

# library load-up ---------------------------------------------------------
library(data.table)
library(InteractionSet)
library(tidyverse)

diff_rnaseq <- readRDS("data/processed/rna/YAPP_HEK_rnaseq_anchors.rds") |> 
  keepStandardChromosomes(pruning.mode = "coarse")

## Reading in loops
## Read in data & format as GInteractions
peaks <- readRDS("data/processed/atac/YAPP_hic_diff_ATACcounts.rds")|> 
  keepStandardChromosomes(pruning.mode = "coarse")

## separate ATAC peaks into gained, lost, static 
gained_peaks <- peaks |>
  subset(padj < 0.1 & log2FoldChange > 0)

lost_peaks <- peaks |>
  subset(padj < 0.1 & log2FoldChange < 0)

static_peaks <- peaks |>
  subset(padj > 0.1)

## Find overlap between unique ATAC anchors & RNA-seq data (using peaks with no normalization)

gained <- subsetByOverlaps(diff_rnaseq, gained_peaks)
lost <- subsetByOverlaps(diff_rnaseq, lost_peaks)
static <- subsetByOverlaps(diff_rnaseq, static_peaks)

## create dataframes of peaks 
gained_df <- as.data.frame(gained)
gained_df$type <- "gained"

lost_df <- as.data.frame(lost)
lost_df$type <- "lost"

static_df <- as.data.frame(static)
static_df$type <- "static"

library(dplyr)

## combine default DESeq2 normalized dataframes
combined1 <- bind_rows(gained_df, lost_df)
combined <- bind_rows(combined1, static_df)

## Plot the RNA log2FoldChange values of "gained", "static", "lost" loops as boxplot
library(ggplot2)

pdf(
  file = "plots/YAPP_HEK_atac_boxplots.pdf",
  height = 3,
  width = 4.5
)
a <- dev.cur()
png(
  file = "vignettes/assets/YAPP_HEK_atac_boxplots.png",
  height = 3,
  width = 4.5
)
dev.control("enable")

ggplot(
  data = combined,
  aes(
    x = reorder(type, desc(type)),
    y = log2FoldChange
  )
) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_boxplot(
    outlier.color = NA,
    fill = "#2E86C1",
    alpha = .5
  ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  ) +
  labs(x = "ATAC regions")

dev.copy(which = a)
dev.off()
dev.off()

## wilcox test for significance
wilcox.test(gained_df$log2FoldChange)

# Create List of genes in gained loops and save as .txt
cat(gained_df$symbol, sep = ",\n")

capture.output(cat(static_df$symbol, sep = ",\n"), file = "data/processed/ATAC/gainedATAC_genes.txt")