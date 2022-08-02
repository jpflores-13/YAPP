
## Gained differential loop visualizations (HiC rectangles)


# Load data ---------------------------------------------------------------
library(plotgardener)
library(hictoolsr)
library(InteractionSet)
library(tidyverse)
library(glue)
library(dbscan)
library(data.table)

load("data/processed/hic/diff_loopCounts_bothDroso.rda")

# Setting gained and lost loops -------------------------------------------

## all loops with a p-value > 0.05
static <- subset(diff_loopCounts, pvalue > 0.05)

## gained loops are loops with a p-value < 0.05 and a (+) log2FC
gained <- subset(diff_loopCounts, pvalue <= 0.05 & log2FoldChange > 0)

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


# create loop lists --------------------------------------------------------
cont_isDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/isDroso/cont"),
                            full.names = T,
                            pattern = "5kbLoops")

cont_noDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/noDroso/cont"),
                            full.names = T,
                            pattern = "5kbLoops")

cont_loops <- c(cont_noDroso_loops, cont_isDroso_loops) |> 
  mergeBedpe(res = 1e4,
             selectCol = 12,
             dist_method = "manhattan",
             minPts = 2)

sorb_isDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/isDroso/sorb",
                                 full.names = T,
                                 pattern = "5kbLoops")

sorb_noDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/noDroso/sorb",
                                 full.names = T,
                                 pattern = "5kbLoops")

sorb_loops <- c(sorb_isDroso_loops, sorb_noDroso_loops) |> 
  mergeBedpe(res = 1e4,
             selectCol = 12,
             dist_method = "manhattan",
             minPts = 2)

omega_isDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/isDroso/omega/",
                                 full.names = T,
                                 pattern = "5kbLoops")

omega_noDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/noDroso/omega/",
                                 full.names = T,
                                 pattern = "5kbLoops")

omega_loops <- c(omega_isDroso_loops, omega_noDroso_loops) |> 
  mergeBedpe(res = 1e4,
             selectCol = 12,
             dist_method = "manhattan",
             minPts = 2)


## use tidyverse to remove `chr` pattern in chromosome{1,2} columns
cont_loops <- 
  cont_loops |> 
  mutate(chromosome1 = str_remove(chromosome1, "chr"),
         chromosome2 = str_remove(chromosome2, "chr"))

sorb_loops <- 
  sorb_loops |> 
  mutate(chromosome1 = str_remove(chromosome1, "chr"),
         chromosome2 = str_remove(chromosome2, "chr"))

omega_loops <- 
  omega_loops |> 
  mutate(chromosome1 = str_remove(chromosome1, "chr"),
         chromosome2 = str_remove(chromosome2, "chr"))

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/gainedLoops_hic_HEK.pdf",
    width = 5.5,
    height = 8)

## Loop through each region
for(i in seq_along(loopRegions_gained)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = gsub('chr', '',loopRegions_gained$seqnames[i]),
                chromstart = loopRegions_gained$start[i],
                chromend = loopRegions_gained$end[i],
                zrange = c(0,100),
                norm = "KR",
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
    plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic", 
                  params = p,
                  y = 0.5) |> 
    annoPixels(data = omega_loops,
               shift = 0.5)
  
  ## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 
  
  control <-
    plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                  params = p,
                  y = 2.75) |> 
    annoPixels(control,
               data = cont_loops,
               shift = 1)
  
  ## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls
  
  sorb <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                  params = p,
                  y = 5) |> 
    annoPixels(sorb, data = sorb_loops, 
               shift = 1)
  
  ## Plot genes
  plotGenes(param = p,
            chrom = paste0('chr', p$chrom),
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
  
  plotText(label = sprintf("Gained YAPP Loop %s", i),
           x = 2.75,
           y = 0.1,
           just = c("center", "top"))
  
}
dev.off()