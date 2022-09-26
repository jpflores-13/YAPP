## How do ATAC peaks change at gained, static, and lost loop anchors? 

# Data import -------------------------------------------------------------

ATAC <- readRDS("data/processed/atac/YAPP_hic_diff_ATACcounts.rds")
loops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

# subset data -------------------------------------------------------------

## separate loops into gained, lost, static
gained <- loops |>
  subset(padj < 0.1 & log2FoldChange > 0)

lost <- loops |>
  subset(padj < 0.1 & log2FoldChange < 0)

static <- loops |>
  subset(padj > 0.1)


# subset by overlap -------------------------------------------------------

## subset gained peaks with gained loops
gained_peak_loops <- subsetByOverlaps(ATAC, gained)
lost_peak_loops <- subsetByOverlaps(ATAC, lost)
static_peak_loops <- subsetByOverlaps(ATAC, static)

## create dataframes of loops
gained_df <- as.data.frame(gained_peak_loops)
gained_df$type <- "gained"


lost_df <- as.data.frame(lost_peak_loops)
lost_df$type <- "lost"

static_df <- as.data.frame(static_peak_loops)
static_df$type <- "static"

library(dplyr)

combined1 <- bind_rows(gained_df, lost_df)
combined <- bind_rows(combined1, static_df)


## Plot the RNA log2FoldChange values of "gained", "static", "lost" loops as boxplot
library(ggplot2)

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
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  ) +
  labs(x = "Loops",
       y = "log2FC ATAC peaks")








