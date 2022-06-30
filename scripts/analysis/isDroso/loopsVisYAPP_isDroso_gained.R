
## Gained differential loop visualization w/ no normalization


# Load data ---------------------------------------------------------------

library(plotgardener)
library(tidyverse)

load("data/output/isDroso/diff_loopCounts.rda")

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

loopRegions_gained <- loopRegions_gained |>
  mutate(seqnames = str_remove(seqnames, "chr"))

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/gainedLoops_isDroso.pdf",
    width = 8.5,
    height = 11)

## Loop through each region
for(i in 1:nrow(loopRegions_gained)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg19",
                resolution = 10e3,
                chrom = loopRegions_gained$seqnames[i],
                chromstart = loopRegions_gained$start[i],
                chromend = loopRegions_gained$end[i],
                zrange = c(0,100),
                norm = "KR",
                x = 0.25,
                width = 3.5,
                length = 3.5,
                height = 3.5,
                fill = "#37a7db",
                linecolor = "#37a7db")


# Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 8.5, height = 5,
             xgrid = 0, ygrid = 0, showGuides = T)
  
  ## Plot top left Hi-C triangle for `-isDroso = TRUE` + annotate
  topLeft <- 
    plotHicSquare(data = "data/hic/dietJuicerMerge/cont/YAPP_HEK_control_inter_30.hic", 
                  params = p,
                  half = "top",
                  y = 0.5) |> 
    annoPixels(data = "data/hic/isDroso_loops/cont/5kbLoops.txt",
               type = "box",
               half = "top", shift = 1)
  
  ## Plot bottom left Hi-C triangle for `-isDroso = TRUE` + annotate
  
  bottomLeft <-
    plotHicSquare(data = "data/hic/dietJuicerMerge/sorb/YAPP_HEK_sorbitol_inter_30.hic",
                  params = p,
                  half = "bottom",
                  y = 0.5) |> 
    annoPixels(data = "data/hic/isDroso_loops/sorb/5kbLoops.txt",
               type = "box",
               half = "bottom",
               shift = 1)
  

  ## Plot top right Hi-C triangle for `-isDroso = FALSE` + annotate

  topRight <- 
    plotHicSquare(data = "data/hic/dietJuicerMerge/cont/YAPP_HEK_control_inter_30.hic", 
                  params = p,
                  half = "top",
                  x = 4.25,
                  y = 0.5) |> 
    annoPixels(data = "data/hic/noDroso_loops/cont/5kbLoops.txt",
               type = "box",
               half = "top", shift = 1)
  
  ## Plot bottom right Hi-C triangle for `-isDroso = FALSE` + annotate
  
  bottomRight <-
    plotHicSquare(data = "data/hic/dietJuicerMerge/sorb/YAPP_HEK_sorbitol_inter_30.hic",
                  params = p,
                  half = "bottom",
                  x = 4.25,
                  y = 0.5) |> 
    annoPixels(data = "data/hic/noDroso_loops/sorb/5kbLoops.txt",
               type = "box",
               half = "bottom",
               shift = 1)
  
  ## Plot genes
  plotGenes(param = p, x = 0.25, y = 4.00, height = 0.5)
  
  ## Plot genes
  plotGenes(param = p, x = 4.25, y = 4.00, height = 0.5)
  
  ## Plot genome label
  plotGenomeLabel(params = p, x = 0.25, y = 4.50)
  
  ## Plot genome label
  plotGenomeLabel(params = p, x = 4.25, y = 4.50)

# Annotate Hi-C triangles by treatment ------------------------------------

  
  ## `-isDroso = TRUE`
  plotText(label = "Untreated",
           x = 0.25,
           y = 0.55,
           just = c("top", "left"))
  
  plotText(label = "Sorbitol",
           x = 3.75,
           y = 3.95,
           just = c("bottom", "right"))
  
  ## `-isDroso = FALSE`
  plotText(label = "Untreated",
           x = 4.25,
           y = 0.55,
           just = c("top", "left"))
  
  plotText(label = "Sorbitol",
           x = 7.75,
           y = 3.95,
           just = c("bottom", "right"))
  

# Annotate Hi-C triangles by `-isDroso` parameter -------------------------

  
  ## annotate columns (untreated)
  plotText(label = "-`isDroso = TRUE`",
           x = 1.5,
           y = 0.25,
           just = c("top", "left"))
  
  ## annotate columns (sorbitol)
  plotText(label = "-`isDroso = FALSE`",
           x = 5.50,
           y = 0.25,
           just = c("top", "left"))
  
  plotText(label = sprintf("Gained YAPP Loop %s", i),
           x = 4.0,
           y = 4.95,
           just = c("center", "bottom"))
  
}

dev.off()

