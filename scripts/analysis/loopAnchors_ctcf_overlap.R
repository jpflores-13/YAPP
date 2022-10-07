## Make pie charts side by side showing % of loop anchors that overlap a CTCF peak 
## (control on left, gained on right)

library(tidyverse)
library(mariner)
library(plyranges)
library(InteractionSet)
library(purrr)
library(plotgardener)
library(RColorBrewer)

## Read in and prune narrowPeak files
peakFiles <- list.files(path = "data/raw/chip/", full.names = TRUE)
chipPeaks <- 
  lapply(peakFiles, read_narrowpeaks) |>
  lapply(keepStandardChromosomes, pruning.mode = 'coarse') 

## Read in loops 
all_loops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

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
                       replace = FALSE)

gained_loops <- all_loops |> 
  subset(padj < 0.1 & log2FoldChange > 0)

chip_gained <- subsetByOverlaps(chipPeaks[[1]], gained_loops)
chip_control <- subsetByOverlaps(chipPeaks[[1]], nullSet)

# data visualization ------------------------------------------------------

# Load ggplot2
library(ggplot2)
library(dplyr)

# Create Data
data <- data.frame(
  group=LETTERS[1:5],
  value=c(13,7,9,21,2)
)

# Compute the position of labels
data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  
  geom_text(aes(y = ypos, label = group), color = "white", size=6) +
  scale_fill_brewer(palette="Set1")
