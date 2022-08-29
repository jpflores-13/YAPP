## ATAC peaks at static, gained, and lost anchors


# library load-up ---------------------------------------------------------
library(data.table)
library(InteractionSet)
library(tidyverse)

# data import -------------------------------------------------------------
atac_peaks <- read.table("data/raw/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv",
                         header = T)

atac_peaks <- GRanges(seqnames = Rle(atac_peaks$chr), 
                      ranges = IRanges(start = atac_peaks$start, end = atac_peaks$stop), 
                      peak = rownames(atac_peaks))

diffLoops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

gained_loops <- diffLoops |> 
  subset(padj < 0.1 & log2FoldChange > 0)

gained_anchors <- subsetByOverlaps(atac_peaks, gained_loops)
  
lost_loops <- diffLoops |> 
  subset(padj < 0.1 & log2FoldChange < 0) 

lost_anchors <- subsetByOverlaps(atac_peaks, lost_loops)

static_loops <- diffLoops |> 
  subset(padj > 0.1)

static_anchors <- subsetByOverlaps(atac_peaks, static_loops)

## create dataframes of peaks
gained_df <- as.data.frame(gained_anchors)
gained_df$type <- "gained"


lost_df <- as.data.frame(lost_anchors)
lost_df$type <- "lost"

static_df <- as.data.frame(static_anchors)
static_df$type <- "static"

## combine data
combined <- bind_rows(gained_df, lost_df, static_df) |> 
  tibble()

combined$peak <- as.integer(combined$peak)
# visualization -----------------------------------------------------------

ggplot(data = combined,
  aes(x = reorder(type, desc(type)),
    y = peak)) +
  geom_hline(yintercept = median(combined$peak), color = "gray", linetype = "dashed") +
  geom_boxplot(outlier.color = NA,
    fill = "#2E86C1",
    alpha = .5) +
  # coord_cartesian(ylim = c(-1, 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "Loop")












