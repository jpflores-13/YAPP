## differential loop visualization w/ no normalization

library(plotgardener)

load("data/output/diff_loopCounts.rda")

print(diff_loopCounts)

## all loops with a p-value > 0.05
static <- subset(diff_loopCounts, pvalue > 0.05)
static

## gained loops are loops with a p-value < 0.05 and a (+) log2FoldChange
gained <- subset(diff_loopCounts, pvalue <= 0.05 & log2FoldChange > 0)
gained

gained_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange > 0)
gained_adj

## filter for the best gained loops
bestGained <- head(diff_loopCounts[order(diff_loopCounts$log2FoldChange, decreasing = T)],50)
bestGained

## lost loops are loops with a p-value < 0.05 and a (-) log2FoldChange
lost <- subset(diff_loopCounts, pvalue <= 0.05 & log2FoldChange < 0)
lost

lost_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange < 0)
lost_adj

## filter for the best lost loops
bestLost <- tail(diff_loopCounts[order(diff_loopCounts$log2FoldChange, decreasing = T)],50)
bestLost

## Create plotting regions from loop anchors
# all gained loops
# loopRegions <- 
#   GRanges(seqnames = as.character(seqnames(anchors(x = gained, "first"))),
#           ranges = IRanges(start = start(anchors(gained, 'first')),
#                            end = end(anchors(gained, 'second'))))

# top 50 gained loops
loopRegions_gained <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = gained_adj, "first"))),
          ranges = IRanges(start = start(anchors(gained_adj, 'first')),
                           end = end(anchors(gained_adj, 'second'))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_gained <- loopRegions_gained + buffer
loopRegions_gained <- as.data.frame(loopRegions_gained)

##make pdf
pdf(file = "plots/old_norm/diffLoopVis_gained.pdf",
    width = 5,
    height = 6)

## Loop through each region
for(i in 1:nrow(loopRegions_gained)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg19",
                resolution = 10e3,
                chrom = loopRegions_gained$seqnames[i],
                chromstart = loopRegions_gained$start[i],
                chromend = loopRegions_gained$end[i],
                zrange = c(0,100),
                norm = "NONE",
                x = 0.5,
                width = 4,
                length = 4,
                height = 2,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  ## Begin visualization
  ## Make page
  pageCreate(width = 5, height = 6,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot Hi-C
  plotHicRectangle(data = "data/omega_hic/GREG_HEK_NUP98PHF23_inter_30.hic",
                   params = p,
                   y = 0.5)
  
  plotHicRectangle(data = "data/omega_hic/LEUK_HEK_PJA27_NAN_NAN_S_0.1.1_megaMap_inter_30.hic",
                   params = p,
                   y = 2.6)
  
  ## Plot Signal track
  plotSignal(data = "data/chip/LEUK_signal/JCH75.bw",
             params = p,
             y = 5,
             height = 0.25,
             just = c("left", "bottom"))
  
  ## Plot genes
  plotGenes(param = p, y = 5.25, height = 0.5)
  
  ## Plot genome label
  plotGenomeLabel(params = p, y = 5.75)
  
  ##annotate
  plotText(label = "PHF23",
           x = 0.55,
           y = 0.55,
           just = c("top", "left"))
  
  plotText(label = "N-IDR_FS/A9",
           x = 0.55,
           y = 2.65,
           just = c("top", "left"))
  
  plotText(label = "CTCF",
           x = 0.55,
           y = 5.05,
           just = c("top", "left"),
           fontsize = 8,
           fontcolor = "#37a7db")
  
  plotText(label = sprintf("PHF23-specific loop %s", i),
           x = 2.5,
           y = 0.375,
           just = c("center", "bottom"))
  
}

dev.off()


########################################################################

# analysis for lost loops

# # all lost loops
# loopRegions <-
#   GRanges(seqnames = as.character(seqnames(anchors(x = lost, "first"))),
#           ranges = IRanges(start = start(anchors(lost, 'first')),
#                            end = end(anchors(lost, 'second'))))

# top 50 gained loops
loopRegions_lost <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = lost_adj, "first"))),
          ranges = IRanges(start = start(anchors(lost_adj, 'first')),
                           end = end(anchors(lost_adj, 'second'))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_lost <- loopRegions_lost + buffer
loopRegions_lost <- as.data.frame(loopRegions_lost)

##make pdf
pdf(file = "plots/old_norm/diffLoopVis_lost.pdf",
    width = 5,
    height = 6)

## Loop through each region
for(i in 1:nrow(loopRegions_lost)){
  
  ## Define parameters
  p <- pgParams(assembly = "hg19",
                resolution = 10e3,
                chrom = loopRegions_lost$seqnames[i],
                chromstart = loopRegions_lost$start[i],
                chromend = loopRegions_lost$end[i],
                zrange = c(0,100),
                norm = "KR",
                x = 0.5,
                width = 4,
                length = 4,
                height = 2,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  ## Begin visualization
  ## Make page
  pageCreate(width = 5, height = 6,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot Hi-C
  plotHicRectangle(data = "data/omega_hic/GREG_HEK_NUP98PHF23_inter_30.hic",
                   params = p,
                   y = 0.5)
  
  plotHicRectangle(data = "data/omega_hic/LEUK_HEK_PJA27_NAN_NAN_S_0.1.1_megaMap_inter_30.hic",
                   params = p,
                   y = 2.6)
  
  ## Plot Signal track
  plotSignal(data = "data/chip/LEUK_signal/JCH75.bw",
             params = p,
             y = 5,
             height = 0.25,
             just = c("left", "bottom"))
  
  ## Plot genes
  plotGenes(param = p, y = 5.25, height = 0.5)
  
  ## Plot genome label
  plotGenomeLabel(params = p, y = 5.75)
  
  ##annotate
  plotText(label = "PHF23",
           x = 0.55,
           y = 0.55,
           just = c("top", "left"))
  
  plotText(label = "N-IDR_FS/A9",
           x = 0.55,
           y = 2.65,
           just = c("top", "left"))
  
  plotText(label = "CTCF",
           x = 0.55,
           y = 5.05,
           just = c("top", "left"),
           fontsize = 8,
           fontcolor = "#37a7db")
  
  plotText(label = sprintf("N-IDR_FS/A9-specific loop %s", i),
           x = 2.5,
           y = 0.375,
           just = c("center", "bottom"))
  
}

dev.off()

