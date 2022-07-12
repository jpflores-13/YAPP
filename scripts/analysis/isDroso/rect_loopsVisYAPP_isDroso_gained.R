
## Gained differential loop visualization w/ no normalization (HiC rectangles)


# Load data ---------------------------------------------------------------

library(plotgardener)
library(data.table)
library(dbscan)
library(hictoolsr)
library(InteractionSet)
library(tidyverse)

load("data/processed/hic/isDroso/diff_loopCounts.rda")

print(diff_loopCounts)

# Setting gained and lost loops -------------------------------------------

## all loops with a p-value > 0.05
static <- subset(diff_loopCounts, pvalue > 0.05)
static

## gained loops are loops with a p-value < 0.05 and a (+) log2FoldChange
gained <- subset(diff_loopCounts, pvalue <= 0.05 & log2FoldChange > 0)
gained

gained_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange > 0)
gained_adj

## filter for the best gained loops
bestGained <- head(gained_adj[order(gained_adj$log2FoldChange, decreasing = T)],100)
bestGained <- head(gained_adj[order(gained_adj$padj, decreasing = F)], 100)
bestGained

## Create plotting regions from loop anchors
# all gained loops
# loopRegions <- 
#   GRanges(seqnames = as.character(seqnames(anchors(x = gained, "first"))),
#           ranges = IRanges(start = start(anchors(gained, 'first')),
#                            end = end(anchors(gained, 'second'))))

# top 50 gained loops
bestGained <- swapAnchors(bestGained)

loopRegions_gained <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
          ranges = IRanges(start = start(anchors(bestGained, "first")),
                           end = end(anchors(bestGained, "second"))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_gained <- loopRegions_gained + buffer
loopRegions_gained <- as.data.frame(loopRegions_gained)

## tidyverse
loopRegions_gained <- loopRegions_gained |>
  mutate(seqnames = str_remove(seqnames, "chr"))

# create omega loop list  -------------------------------------------------

omega_loops <- c("data/raw/hic/hg38/sip-loops/isDroso/omega/5kbLoops.txt",
                "data/raw/hic/hg38/sip-loops/isDroso/sorb/5kbLoops.txt",
                "data/raw/hic/hg38/sip-loops/isDroso/cont/5kbLoops.txt")

# omega_loops <-
#   omega_loops |>
#   lapply(fread) |>
#   {\(x) do.call(rbind, x)}()


# get rid of 'chrom' ------------------------------------------------------

# omega_loops <- omega_loops |> 
#   mutate(chromosome1 = str_remove(chromosome1, "chr"))
# 
# omega_loops <- omega_loops |> 
#   mutate(chromosome2 = str_remove(chromosome2, "chr"))

# subset loop files -------------------------------------------------------
# cont_loops <- "data/hic/isDroso/cont/loops/5kbLoops.txt"
# 
# cont_loops <-
#   cont_loops |>
#   lapply(fread) |>
#   {\(x) do.call(rbind, x)}()
# 
# cont_loops <- cont_loops |> 
#   mutate(chromosome1 = str_remove(chromosome1, "chr"))
# 
# cont_loops <- cont_loops |> 
#   mutate(chromosome2 = str_remove(chromosome2, "chr"))
# 
# sorb_loops <- "data/hic/isDroso/sorb/loops/5kbLoops.txt"
# sorb_loops <-
#   sorb_loops |>
#   lapply(fread) |>
#   {\(x) do.call(rbind, x)}()
# 
# sorb_loops <- sorb_loops |> 
#   mutate(chromosome1 = str_remove(chromosome1, "chr"))
# 
# sorb_loops <- sorb_loops |> 
#   mutate(chromosome2 = str_remove(chromosome2, "chr"))

# get rid of 'chrom' for noDroso files ------------------------------------

# cont_loops_noDroso <- "data/hic/noDroso/cont/loops/5kbLoops.txt"
# cont_loops_noDroso <-
#   cont_loops_noDroso |>
#   lapply(fread) |>
#   {\(x) do.call(rbind, x)}()
# 
# cont_loops_noDroso <- cont_loops_noDroso |> 
#   mutate(chromosome1 = str_remove(chromosome1, "chr"))
# 
# cont_loops_noDroso <- cont_loops_noDroso |> 
#   mutate(chromosome2 = str_remove(chromosome2, "chr"))
# 
# sorb_loops_noDroso <- "data/hic/noDroso/sorb/loops/5kbLoops.txt"
# sorb_loops_noDroso <-
#   sorb_loops_noDroso |>
#   lapply(fread) |>
#   {\(x) do.call(rbind, x)}()
# 
# sorb_loops_noDroso <- sorb_loops_noDroso |> 
#   mutate(chromosome1 = str_remove(chromosome1, "chr"))
# 
# sorb_loops_noDroso <- sorb_loops_noDroso |> 
#   mutate(chromosome2 = str_remove(chromosome2, "chr"))


# tmp ---------------------------------------------------------------------


# Subset Loop Overlaps ----------------------------------------------------

# ## control loops
# loopOverlaps_cont <- c("data/hic/isDroso/cont/loops/5kbLoops.txt",
#                           "data/hic/noDroso/cont/loops/5kbLoops.txt")
# ## sorbitol loops
# loopOverlap_sorb <- c("data/hic/isDroso/sorb/loops/5kbLoops.txt",
#                        "data/hic/noDroso/sorb/loops/5kbLoops.txt")
# 
# ## convert to GI
# loopOverlap_merge <- 
#   mergeBedpe(bedpeFiles = c(loopOverlaps_cont,loopOverlap_sorb),
#                                 res = 10000,
#                                 selectCol = 12,
#                                 dist_method = "manhattan",
#                                 minPts = 2) |> 
#   as_ginteractions() |> 
#   binBedpe(res = 10e3, 
#            a1Pos = "center", 
#            a2Pos = "center") 

## subset by overlaps

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/gainedLoops_isDroso_rect.pdf",
    width = 7,
    height = 8)

## Loop through each region
for(i in 1:nrow(loopRegions_gained)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg19",
                resolution = 5e3,
                chrom = loopRegions_gained$seqnames[1],
                chromstart = loopRegions_gained$start[1],
                chromend = loopRegions_gained$end[1],
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
  pageCreate(width = 7, height = 8,
             xgrid = 0, ygrid = 0, showGuides = T)
  
  ## Plot top omega Hi-C rectangle + annotate
  top <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic", 
                  params = p,
                  y = 0.5)
  
  annoPixels(plot = top,
             data = "data/raw/hic/hg38/sip-loops/isDroso/omega/5kbLoops.txt",
             type = "box",
             shift = 1)
  
  ## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 
  
  middle <-
    plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                  params = p,
                  y = 2.75) 
  
    annoPixels(plot = middle,
               data = cont_loops_noDroso,
               type = "box",
               shift = 1)
  
    annoPixels(plot = middle,
               data = cont_loops,
               type = "circle",
               shift = 1)
    
  ## for overlap boxes, read data in, convert to GInteractions, subset by overlaps
  
  
  ## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls
  
  bottom <- 
    plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                  params = p,
                  y = 5)
  
    annoPixels(plot = bottom,
             data = sorb_loops_noDroso,
             type = "box",
             shift = 1)
    
    annoPixels(plot = bottom,
               data = sorb_loops,
               type = "circle",
               shift = 1)

  ## Plot genes
  plotGenes(param = p, x = 0.25, y = 7.25, height = 0.5)
  
  ## Plot genome label
  plotGenomeLabel(params = p, x = 0.25, y = 7.75)
  
  # Annotate Hi-C rectangles by treatment ------------------------------------

  plotText(label = "omega",
           x = 0.25,
           y = 0.35,
           just = c("top", "left"))
  
  plotText(label = "untreated",
           x = 0.25,
           y = 2.6,
           just = c("top", "left"))
  
  plotText(label = "sorbitol",
           x = 0.25,
           y = 4.85,
           just = c("top", "left"))
  
  plotText(label = "`-isDroso = FALSE`",
           y = 3.5,
           x = 6,
           just = c("center", "center"),
           col = "red")
  
  plotText(label = "`-isDroso = TRUE`",
           y = 4,
           x = 6,
           just = c("center", "center"),
           col = "green")
  
  plotText(label = sprintf("Gained YAPP Loop %s", i),
           x = 2.75,
           y = 0.1,
           just = c("center", "top"))
  
}

dev.off()

