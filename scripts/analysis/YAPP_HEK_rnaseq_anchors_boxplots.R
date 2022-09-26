library(InteractionSet)

diff_rnaseq <- readRDS("data/processed/rna/YAPP_HEK_rnaseq_anchors.rds")

## Reading in loops
## Read in data & format as GInteractions

loops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

## separate loops into gained, lost, static
gained <- loops |>
  subset(padj < 0.1 & log2FoldChange > 0)

lost <- loops |>
  subset(padj < 0.1 & log2FoldChange < 0)

static <- loops |>
  subset(padj > 0.1)


## Find overlap between unique loop anchors

gained <- subsetByOverlaps(diff_rnaseq, gained)
lost <- subsetByOverlaps(diff_rnaseq, lost)
static <- subsetByOverlaps(diff_rnaseq, static)


## create dataframes of loops
gained_df <- as.data.frame(gained)
gained_df$type <- "gained"


lost_df <- as.data.frame(lost)
lost_df$type <- "lost"

static_df <- as.data.frame(static)
static_df$type <- "static"

library(dplyr)

combined1 <- bind_rows(gained_df, lost_df)
combined <- bind_rows(combined1, static_df)


## Plot the RNA log2FoldChange values of "gained", "static", "lost" loops as boxplot
library(ggplot2)

pdf(
  file = "plots/YAPP_HEK_rnaseq_anchors_boxplots.pdf",
  height = 3,
  width = 4.5
)
a <- dev.cur()
png(
  file = "vignettes/assets/YAPP_HEK_rnaseq_anchors_boxplots.png",
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
  labs(x = "Loops")

dev.copy(which = a)
dev.off()
dev.off()

## wilcox test for significance
wilcox.test(gained_df$log2FoldChange)

# Create List of genes in gained loops and save as .txt
cat(gained_df$symbol, sep = ",\n")

capture.output(cat(static_df$symbol, sep = ",\n"), file = "data/processed/rna/gainedLoop_genes.txt")
