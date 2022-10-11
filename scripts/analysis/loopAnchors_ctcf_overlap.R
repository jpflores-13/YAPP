## Make pie charts side by side showing % of loop anchors that overlap a CTCF peak 
## (control on left, gained on right)

library(tidyverse)
library(mariner)
library(plyranges)
library(InteractionSet)
library(purrr)
library(plotgardener)
library(RColorBrewer)
library(nullranges)

## Read in and prune narrowPeak files
peakFiles <- list.files(path = "data/raw/chip/", full.names = TRUE)
chipPeaks <- 
  map(peakFiles, read_narrowpeaks) |>
  map(keepStandardChromosomes, pruning.mode = 'coarse')

## Read in loops 
all_loops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds") |> 
  reduceRegions() |> 
  regions()
  
## create metadata columns for contact frequency & size 
mcols(all_loops)$loop_size <- pairdist(all_loops)
mcols(all_loops)$loop_type <- case_when(
  mcols(all_loops)$padj < 0.1 & mcols(all_loops)$log2FoldChange > 0 ~ "gained",
  mcols(all_loops)$padj < 0.1 & mcols(all_loops)$log2FoldChange < 0 ~ "lost",
  mcols(all_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

## contact frequency
mcols(all_loops)$sorb_contacts <- mcols(all_loops)$YAPP_HEK_sorbitol_4_2_inter_30.hic +
  mcols(all_loops)$YAPP_HEK_sorbitol_5_2_inter_30.hic + mcols(all_loops)$YAPP_HEK_sorbitol_6_2_inter_30.hic

## use matchRanges to select a null set of control sample loops that is matched for size & contact frequency
nullSet <- matchRanges(focal = all_loops[mcols(all_loops)$loop_type == "gained"],
                       pool = all_loops[!mcols(all_loops)$loop_type == "gained"],
                       covar = ~ loop_size + sorb_contacts, 
                       method = 'stratified',
                       replace = FALSE) |> 
  reduceRegions() |> 
  regions()

gained_loops <- all_loops |> 
  subset(padj < 0.1 & log2FoldChange > 0) |> 
  reduceRegions() |> 
  regions()

## returns only the ranges in the first object that have overlaps with any ranges in the second object.

chip_gained <- countOverlaps(gained_loops, chipPeaks[[1]]) 

chip_nullSet <- countOverlaps(nullSet, chipPeaks[[1]]) 

chip_control <- countOverlaps(chipPeaks[[1]], all_loops)

## Calculate percentages of the total number of loops and round


# data visualization ------------------------------------------------------

ggplot(data, aes(x="", y=amount, fill=category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(amount, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_brewer(palette="Blues")







###############notes
bothP       <- round(bothN/totalN*100)
oneP        <- round(oneN/totalN*100)
neitherP    <- round(neitherN/totalN*100)


## Create data frame to plot results
df  <- data.frame(
  group = c("Both", "One", "Neither"),
  number = c(bothN, oneN, neitherN),
  percent = c(bothP, oneP, neitherP)
)

## Plot results
library(ggplot2)
ggplot(data = df, aes(x = 1, y = percent, fill = group))+
  geom_col(col = "white")+
  coord_polar("y") +
  scale_fill_manual(values = c("#2171b5", "#bdd7e7", "#6baed6"))+
  labs(title = "Loop anchors with bound CTCF")+
  annotate(geom = "text",
           x = c(1.0, 1.25, 1.75),
           y = c(60, 8, 17.25),
           label = paste0(df$group, " (", df$percent, "%)"),
           col = c("white", "white", "black"),
           size = 4.5)+
  theme_void()+
  theme(
    plot.title = element_text(hjust = 0.5, vjust = -15, face = "bold"),
    legend.position = "none", text = element_text(face = "bold")
  )