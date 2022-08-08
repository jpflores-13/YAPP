## Gained Loop APA plots for HEK Hi-C

library(tidyverse)
library(strawr)
library(hictoolsr)
library(GenomeInfoDb)
library(RColorBrewer)
library(plotgardener)

# Load in all loops -------------------------------------------------------
diff_loopCounts <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

## convert seqlevelStyles
seqlevelsStyle(diff_loopCounts) <- "UCSC"

## Filter all static loops
staticLoops <- diff_loopCounts |> 
  subset(padj > 0.1)

## Filter out gained loops
gainedLoops <- diff_loopCounts |> 
  subset(padj < 0.1 & log2FoldChange > 0)

## Filter out lost loops
lostLoops <- diff_loopCounts |> 
  subset(padj < 0.1 & log2FoldChange < 0)

## Compile loop lists
loopList <- 
  list(staticLoops = staticLoops,
       gainedLoops = gainedLoops,
       lostLoops = lostLoops)

## Define resolution, buffer, and normalization (pixels from center)
res <- 10e3
buffer <- 10
norm <- "SCALE"

## Filter out short loop interactions
filteredLoops <- 
  map(loopList,filterBedpe, res = res, buffer = buffer) |> 
  'names<-'(value = names(loopList))

## Summary loops
map(filteredLoops, summary)

## Micro-C file paths 
cont_micro_YAPP <- "data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic"
sorb_microc_YAPP <- "data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic"

## Calculate APA matrices for loops from control HEK Micro-C data
cont_APA_mat <- map(filteredLoops, calcApa, hic = cont_micro_YAPP, norm = norm, buffer = buffer)

## Calculate APA matrices for loops from sorbitol HEK Micro-C data
sorb_APA_mat <- map(filteredLoops, calcApa, hic = sorb_microc_YAPP, norm = norm, buffer = buffer)

## Get the number of loops for each condition
nLoops <- map(filteredLoops, length)

## Divide each matrix by nLoops
cont_APA_mat <- Map("/", cont_APA_mat, nLoops)
sorb_APA_mat <- Map("/", sorb_APA_mat, nLoops)

## combine APA matrices to pull out the max value for zrange max
mats_combined <- c(cont_APA_mat,
                   sorb_APA_mat)


# Create plotgardener page ------------------------------------------------
pdf(file = "plots/YAPP_HEK_microc_APA.pdf",
    height = 3,
    width = 4.5)

## Initiate plotgardener page
pageCreate(width = 4.5, height = 3, showGuides = F)

## Define shared parameters
p <- pgParams(x = 0.5,
              y = 0.5,
              width = 1,
              height = 1,
              space = 0.075,
              zrange = c(0, max(unlist(mats_combined))),
              palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))

## Plot gained control APAs
leftTop <- plotApa(params = p, cont_APA_mat$gainedLoops) |> 
  annoHeatmapLegend(x = 1.55, y = 0.5, width = 0.1, 
                    height = 0.5)

## Plot gained treated APAs
leftBottom <- plotApa(params = p, sorb_APA_mat$gainedLoops, y = 1.55) |> 
  annoHeatmapLegend(x = 1.55, y = 1.55, width = 0.1, height = 0.5)

## Plot lost control APAs
middleTop <- plotApa(params = p, cont_APA_mat$lostLoops, x = 1.75) |> 
  annoHeatmapLegend(x = 2.85, y = 0.5, width = 0.1, 
                    height = 0.5)

## Plot lost treated APAs
middleBottom <- plotApa(params = p, sorb_APA_mat$lostLoops, x = 1.75, y = 1.55) |> 
  annoHeatmapLegend(x = 2.85, y = 1.55, width = 0.1, 
                    height = 0.5)  

## Plot static control APAs
topRight <- plotApa(params = p, cont_APA_mat$staticLoops, x = 3.05, y = 0.5) |> 
  annoHeatmapLegend(x = 4.1, y = 0.5, width = 0.1, 
                    height = 0.5)

## Plot static treated APAs
bottomRight <- plotApa(params = p, sorb_APA_mat$staticLoops, x = 3.05, y = 1.55) |> 
  annoHeatmapLegend(x = 4.1, y = 1.55, width = 0.1, height = 0.5)

## Add text labels
plotText(label = "Gained",
         x = 1,
         y = 0.25,
         just = c('center', 'top'))

plotText(label = "Lost",
         x = 2.25,
         y = 0.25,
         just = c('center', 'top'))

plotText(label = "Static",
         x = 3.5,
         y = 0.25,
         just = c('center', 'top'))

plotText(label = "Cont",
         x = 0.45,
         y = 1,
         rot = 90,
         just = c('center', 'bottom'))  

plotText(label = "Sorb",
         x = 0.45,
         y = 2,
         rot = 90,
         just = c('center', 'bottom'))  


dev.off()
