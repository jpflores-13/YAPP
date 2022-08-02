
## Lost differential loop visualization w/ no normalization (HiC Squares/Tris)


# Load Data ---------------------------------------------------------------

library(plotgardener)

load("data/output/isDroso/diff_loopCounts.rda")

print(diff_loopCounts)


# Setting gained and lost loops -------------------------------------------

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
bestLost <- head(lost_adj[order(lost_adj$padj, decreasing = F)], 100)
bestLost

# # all lost loops
# loopRegions <-
#   GRanges(seqnames = as.character(seqnames(anchors(x = lost, "first"))),
#           ranges = IRanges(start = start(anchors(lost, 'first')),
#                            end = end(anchors(lost, 'second'))))

# top 50 lost loops
bestLost <- swapAnchors(bestLost)

loopRegions_lost <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestLost, "first"))),
          ranges = IRanges(start = start(anchors(bestLost, 'first')),
                           end = end(anchors(bestLost, 'second'))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_lost <- loopRegions_lost + buffer
loopRegions_lost <- as.data.frame(loopRegions_lost)

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/lostLoops_isDroso.pdf",
    width = 8.5,
    height = 11)

## Loop through each region
for(i in 1:nrow(loopRegions_lost)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg19",
                resolution = 5e3,
                chrom = loopRegions_lost$seqnames[i],
                chromstart = loopRegions_lost$start[i],
                chromend = loopRegions_lost$end[i],
                zrange = c(0,100),
                norm = "KR",
                x = 0.5,
                width = 3.5,
                length = 3.5,
                height = 3.5,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 8.5, height = 5,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot top Hi-C triangle for `-isDroso = TRUE` + annotate
  topTop <- 
    plotHicSquare(data = "data/hic/isDroso/cont/YAPP_HEK_control_inter_30.hic", 
                  params = p,
                  half = "top",
                  y = 0.5) |> 
    annoPixels(data = "data/hic/isDroso/cont/loops/5kbLoops.txt",
               type = "box",
               half = "top", shift = 1)
  
  ## Plot bottom Hi-C triangle for `-isDroso = TRUE` + annotate
  
  bottomTop <-
    plotHicSquare(data = "data/hic/isDroso/sorb/YAPP_HEK_sorbitol_inter_30.hic",
                  params = p,
                  half = "bottom",
                  y = 0.5) |> 
    annoPixels(data = "data/hic/isDroso/sorb/loops/5kbLoops.txt",
               type = "box",
               half = "bottom",
               shift = 1)
  
  ## Plot top Hi-C triangle for `-isDroso = FALSE` + annotate
  
  topTop_false <- 
    plotHicSquare(data = "data/hic/noDroso/cont/YAPP_HEK_control_inter_30.hic", 
                  params = p,
                  half = "top",
                  x = 4.25,
                  y = 0.5) |> 
    annoPixels(data = "data/hic/noDroso/cont/loops/5kbLoops.txt",
               type = "box",
               half = "top", shift = 1)
  
  ## Plot bottom Hi-C triangle for `-isDroso = FALSE` + annotate
  
  bottomTop_false <-
    plotHicSquare(data = "data/hic/noDroso/sorb/YAPP_HEK_sorbitol_inter_30.hic",
                  params = p,
                  half = "bottom",
                  x = 4.25,
                  y = 0.5) |> 
    annoPixels(data = "data/hic/noDroso/sorb/loops/5kbLoops.txt",
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
  
  plotText(label = sprintf("Lost YAPP Loop %s", i),
           x = 4.0,
           y = 4.95,
           just = c("center", "bottom"))
  
}

dev.off()