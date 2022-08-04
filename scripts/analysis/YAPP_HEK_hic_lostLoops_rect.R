
## Lost differential loop visualizations (HiC rectangles)

# Load data ---------------------------------------------------------------
library(plotgardener)
library(hictoolsr)
library(InteractionSet)
library(tidyverse)
library(glue)
library(dbscan)
library(data.table)

diff_loopCounts <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

# Setting static and lost loops -------------------------------------------

## all loops with a p-value > 0.05
static <- subset(diff_loopCounts, pvalue > 0.05)
static

## lost loops are loops with a p-value < 0.05 and a (-) log2FoldChange
lost <- subset(diff_loopCounts, pvalue <= 0.05 & log2FoldChange < 0)
lost

lost_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange < 0)
lost_adj

## filter for the best lost loops
bestLost <- head(lost_adj[order(lost_adj$log2FoldChange, decreasing = T)],100)
bestLost <- head(lost_adj[order(lost_adj$pvalue, decreasing = F)], 100)
bestLost

# # all lost loops
# loopRegions <-
#   GRanges(seqnames = as.character(seqnames(anchors(x = lost, "first"))),
#           ranges = IRanges(start = start(anchors(lost, 'first')),
#                            end = end(anchors(lost, 'second'))))

# top 50 lost loops
# bestLost <- swapAnchors(bestLost)

loopRegions_lost <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestLost, "first"))),
          ranges = IRanges(start = start(anchors(bestLost, 'first')),
                           end = end(anchors(bestLost, 'second'))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_lost <- loopRegions_lost + buffer
loopRegions_lost <- as.data.frame(loopRegions_lost)

## Use tidyverse to remove `chr`in seqnames
loopRegions_lost <- loopRegions_lost |>
  mutate(seqnames = str_remove(seqnames, "chr"))

# load loop lists --------------------------------------------------------

cont_loops <- readRDS("data/processed/hic/cont_bothDroso_loops.rds")

sorb_loops <- readRDS("data/processed/hic/sorb_bothDroso_loops.rds")

omega_loops <- readRDS("data/processed/hic/omega_bothDroso_loops.rds")

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/YAPP_HEK_hic_lostLoops_rect.pdf",
    width = 5.5,
    height = 8)

## Loop through each region
for(i in 1:nrow(loopRegions_lost)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = paste0("chr",loopRegions_lost$seqnames[i]),
                chromstart = loopRegions_lost$start[i],
                chromend = loopRegions_lost$end[i],
                zrange = c(0,100),
                norm = "SCALE",
                x = 0.25,
                width = 5,
                length = 5,
                height = 2,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 5.5, height = 8,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot top omega Hi-C rectangle + annotate
  omega <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic", 
                     params = p,
                     y = 0.5) |> 
    annoPixels(data = omega_loops,
               shift = 0.5)
  
  ## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 
  
  control <-
    plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                     params = p,
                     y = 2.75) |> 
    annoPixels(data = cont_loops,
               shift = 1)
  
  ## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls
  
  sorb <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                     params = p,
                     y = 5) |> 
    annoPixels(data = sorb_loops, 
               shift = 1)
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 7.25,
            height = 0.5)
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 7.75)
  
  # Annotate Hi-C rectangles by treatment ------------------------------------
  
  plotText(label = "omega",
           x = 0.25,
           y = 0.35,
           just = c("top", "left"))
  
  plotText(label = "control",
           x = 0.25,
           y = 2.6,
           just = c("top", "left"))
  
  plotText(label = "sorbitol",
           x = 0.25,
           y = 4.85,
           just = c("top", "left"))
  
  plotText(label = sprintf("Lost YAPP Loop %s", i),
           x = 2.75,
           y = 0.1,
           just = c("center", "top"))
  
}
dev.off()