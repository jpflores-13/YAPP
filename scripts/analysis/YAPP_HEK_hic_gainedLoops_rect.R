
## Gained differential loop visualizations (HiC rectangles)


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

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange > 0)
gained_adj

## filter for the 100 best gained loops
bestGained <- head(gained_adj[order(gained_adj$log2FoldChange, decreasing = T)],100)
bestGained <- head(gained_adj[order(gained_adj$padj, decreasing = F)], 100)
bestGained

# top 100 gained loops
## If the below function doesn't work, might need to use swapAchors()

loopRegions_gained <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
          ranges = IRanges(start = start(anchors(bestGained, "first")),
                           end = end(anchors(bestGained, "second"))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_gained <- loopRegions_gained + buffer
loopRegions_gained <- as.data.frame(loopRegions_gained)

## Use tidyverse to remove `chr`in seqnames
loopRegions_gained <- loopRegions_gained |>
  mutate(seqnames = str_remove(seqnames, "chr"))

# load loop lists --------------------------------------------------------

cont_loops <- readRDS("data/processed/hic/cont_bothDroso_loops.rds")

sorb_loops <- readRDS("data/processed/hic/sorb_bothDroso_loops.rds")

omega_loops <- readRDS("data/processed/hic/omega_bothDroso_loops.rds")

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/YAPP_HEK_hic_gainedLoops_rect.pdf",
    width = 7,
    height = 8)

## Loop through each region
for(i in 1:nrow(loopRegions_gained)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = paste0("chr",loopRegions_gained$seqnames[i]),
                chromstart = loopRegions_gained$start[i],
                chromend = loopRegions_gained$end[i],
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
  pageCreate(width = 7, height = 8,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot top omega Hi-C rectangle + annotate
  omega <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic", 
                  params = p,
                  y = 0.5)
  ## annotate all pixels called with both `-isDroso true` and `-isDroso false` parameters
    annoPixels(omega, data = omega_loops,
               type = "arrow",
               shift = 0.5)
    
  ## Plot middle Hi-C rectangle 
  
  control <-
    plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                  params = p,
                  y = 2.75)
  
  ## annotate all pixels called with both `-isDroso true` and `-isDroso false` parameters
  annoPixels(control, data = omega_loops,
             shift = 0.5, 
             type = "arrow",
             col = "black")
  
  ## annotate all pixels called with `-isDroso true`
  annoPixels(control, data = "data/raw/hic/hg38/sip-loops/isDroso/cont/5kbLoops.txt",
             shift = 0.5,
             type = "arrow",
             col = "green")
  
  ## annotate all pixels called with `-isDroso false`
  annoPixels(control, data = "data/raw/hic/hg38/sip-loops/noDroso/cont/5kbLoops.txt",
             shift = 0.5, 
             type = "arrow",
             col = "red")
    
  ## Plot bottom Hi-C rectangle 
  
  sorb <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                  params = p,
                  y = 5)
  
  ## annotate all pixels called with both `-isDroso true` and `-isDroso false` parameters
    annoPixels(sorb, data = omega_loops, 
               shift = 0.5)
    
    ## annotate all pixels called with `-isDroso true`
    annoPixels(sorb, data = "data/raw/hic/hg38/sip-loops/isDroso/sorb/5kbLoops.txt",
               shift = 0.5, 
               type = "arrow",
               col = "green")
    
    ## annotate all pixels called with `-isDroso false`
    annoPixels(sorb, data = "data/raw/hic/hg38/sip-loops/noDroso/sorb/5kbLoops.txt",
               shift = 0.5, 
               type = "arrow",
               col = "red")
  
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
  
  plotLegend(legend = c("-isDroso true", "-isDroso false", "omega loops"),
    fill = c("green", "red", "black"),
    border = FALSE,
    x = 5.5, y = 3, width = 1.5, height = 0.7,
    just = c("left", "top"),
    default.units = "inches")
  
  plotText(label = sprintf("Gained YAPP Loop %s", i),
           x = 2.75,
           y = 0.1,
           just = c("center", "top"))
  
}
dev.off()