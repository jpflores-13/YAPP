## YAPP (HEK) Micro-C data by APA plots 

library(tidyverse)
library(strawr)
library(hictoolsr)
library(GenomeInfoDb)
library(RColorBrewer)
library(plotgardener)

## Load in all loops
load("data/processed/microc/YAPP_microc_loopCounts.rda")

contLoops <- readRDS("data/processed/microc/cont_bothDroso_loops.rds") |> 
  as_ginteractions() 

sorbLoops <- readRDS("data/processed/microc/sorb_bothDroso_loops.rds") |> 
  as_ginteractions()

seqlevelsStyle(loopCounts) <- "UCSC"
seqlevelsStyle(contLoops) <- "UCSC"
seqlevelsStyle(sorbLoops) <- "UCSC"

## Compile all YAPP loops and plop them into a list
loopList <- 
  list(allLoops = loopCounts,
       contLoops = contLoops,
       sorbLoops = sorbLoops)

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

## Hi-C file paths 
cont_hic <- "data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic"
sorb_hic <- "data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic"

## Calculate APA matrices for loops from control Hi-C data
cont_APA_mat <- 
  map(filteredLoops, calcApa, hic = cont_hic, norm = norm, buffer = buffer)

## Calculate APA matrices for loops from sorbitol Hi-C data
sorb_APA_mat <-
  map(filteredLoops, calcApa, hic = sorb_hic, norm = norm, buffer = buffer)

## had to debug `calcAPA()` with Eric. Changes made to JP's forked version of hictoolsr
# debugonce(calcApa)

## Get the number of loops for each condition
nLoops <- map(filteredLoops, length)

## Divide each matrix by nLoops
loopApaContHic <- Map("/", cont_APA_mat, nLoops)
loopApaSorbHic <- Map("/", sorb_APA_mat, nLoops)

## plot
plotApa(apa = loopApaContHic$allLoops,
        palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))


# Create plotgardener page ------------------------------------------------
pdf(file = "plots/YAPP_HEK_microc_APA.pdf",
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
              zrange = c(0, max(unlist(c(loopApaContHic, loopApaSorbHic)))),
              palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))

## Define grid of coordinate positions
xpos <- c(p$x, p$x + p$width + p$space, p$x + (p$width + p$space)*2)
ypos <- c(p$y, p$y + p$height + p$space, p$y + (p$height + p$space)*2)

## Plot row of cont APAs
cont_plots <- 
  pmap(list(loopApaContHic, xpos, ypos[1]), \(a, x, y) {
    plotApa(params = p, apa = a, x = x, y = y)
  })

## Plot row of sorbitol APAs
sorb_plots <- 
  pmap(list(loopApaSorbHic, xpos, ypos[2]), \(a, x, y) {
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
plotText(label = c("All loops", "Control loops", "Sorbitol loops"),
         x = xpos + p$width / 2,
         y = ypos[1] - p$space,
         just = c('center', 'bottom'))

plotText(label = c("Cont", "Sorb"),
         x = xpos[1] - p$space,
         y = ypos[1:2] + p$height / 2,
         rot = 90,
         just = c('center', 'bottom'))  
dev.off()






