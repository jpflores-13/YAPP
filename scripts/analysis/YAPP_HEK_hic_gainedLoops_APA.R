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

## Filter out gained control loops
gainedLoops <- diff_loopCounts |> 
  subset(padj < 0.1 & log2FoldChange > 0)

## Define resolution, buffer, and normalization (pixels from center)
res <- 10e3
buffer <- 10
norm <- "SCALE"

## Filter out short loop interactions
filteredLoops <- gainedLoops |> 
  filterBedpe(res = res, 
              buffer = buffer) |> 
  'names<-'(value = names(gainedLoops))

## Summary loops
filteredLoops |> 
  summary()

## Hi-C file paths 
cont_hic_YAPP <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic"
sorb_hic_YAPP <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic"

## Calculate APA matrices for loops from control HEK Hi-C data
cont_APA_mat_hic_YAPP <- filteredLoops |> 
  calcApa(hic = cont_hic_YAPP, norm = norm, buffer = buffer)

## Calculate APA matrices for loops from sorbitol HEK Hi-C data
sorb_APA_mat_hic_YAPP <- filteredLoops |> 
  calcApa(hic = sorb_hic_YAPP, norm = norm, buffer = buffer)

## had to debug `calcAPA()` with Eric. Changes made to JP's forked version of hictoolsr
# debugonce(calcApa)

## Get the number of loops for each condition
nLoops <- filteredLoops |> 
  length()

## Divide each matrix by nLoops
cont_APA_mat_hic_YAPP <- (cont_APA_mat_hic_YAPP/nLoops)
sorb_APA_mat_hic_YAPP <- (sorb_APA_mat_hic_YAPP/nLoops)

## combine APA matrices to pull out the max value for zrange max
mats_combined <- c(cont_APA_mat_hic_YAPP,
                   sorb_APA_mat_hic_YAPP)


# Create plotgardener page ------------------------------------------------
pdf(file = "plots/YAPP_HEK_hic_gainedLoops_APA.pdf",
    height = 3,
    width = 4.25)

## Initiate plotgardener page
pageCreate(width = 4.25, height = 3, showGuides = F)

## Define shared parameters
p <- pgParams(x = 0.5,
              y = 0.5,
              width = 1,
              height = 1,
              space = 0.075,
              zrange = c(0, max(mats_combined)),
              palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))


## Define grid of coordinate positions
xpos <- c(p$x, p$x + p$width + p$space, p$x + (p$width + p$space)*2)
ypos <- c(p$y, p$y + p$height + p$space, p$y + (p$height + p$space)*2)

## Plot row of cont APAs
cont_plot <-
  pmap(list(cont_APA_mat_hic_YAPP, xpos, ypos[1]), \(a, x, y) {
    plotApa(params = p, apa = a, x = x, y = y)
  })

## Plot row of sorbitol APAs
treated_plot <-
  pmap(list(sorb_APA_mat_hic_YAPP, xpos, ypos[2]), \(a, x, y) {
    plotApa(params = p, apa = a, x = x, y = y)
  })

## Add legend
annoHeatmapLegend(plot = cont_plots[[1]],
                  x = p$x + (p$width + p$space)*3,
                  y = ypos[1],
                  width = p$space,
                  height = p$height*0.75,
                  fontcolor = 'black')

## Add text labels
plotText(label = c("HEK Hi-C", "HEK Micro-C", "T47D Hi-C"),
         x = xpos + p$width / 2,
         y = ypos[1] - p$space,
         just = c('center', 'bottom'))

plotText(label = c("Cont", "Treated"),
         x = xpos[1] - p$space,
         y = ypos[1:2] + p$height / 2,
         rot = 90,
         just = c('center', 'bottom'))  
dev.off()
