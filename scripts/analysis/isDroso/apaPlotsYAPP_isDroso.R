## install
# remotes::install_github("EricSDavis/hictoolsr")

## load packages
library(hictoolsr)
library(data.table)
library(tidyverse)
library(InteractionSet)
library(dbscan)

# load data
load("data/output/diff_loopCounts.rda")

# define parameters
res <- 10e3
buffer <- 10
norm <- "KR"
pvthresh <- 0.1

GREG_hic <- "data/omega_hic/GREG_HEK_NUP98PHF23_inter_30.hic"
LEUK_hic <- "data/omega_hic/LEUK_HEK_PJA27_NAN_NAN_S_0.1.1_megaMap_inter_30.hic"


# filter out short loops
diff_loopCounts <- filterBedpe(bedpe = diff_loopCounts,
                               res = res,
                               buffer = buffer)

# pull out what you need 
## all loops with a p-value > 0.05
static <- subset(diff_loopCounts, padj > pvthresh)

## gained loops are loops with a p-value < pvthresh and a (+) log2FoldChange
gained <- subset(diff_loopCounts, padj <= pvthresh & log2FoldChange > 0)

## lost loops are loops with a p-value < pvthresh and a (-) log2FoldChange
lost <- subset(diff_loopCounts, padj <= pvthresh & log2FoldChange < 0)

# calculate APA using .hic files imported from LongLeaf
staticGregAPA <- calcApa(bedpe = static, 
                         hic = GREG_hic,
                         norm = norm,
                         res = res,
                         buffer = buffer) 

staticLeukAPA <- calcApa(bedpe = static, 
                         hic = LEUK_hic,
                         norm = norm,
                         res = res,
                         buffer = buffer)

gainedGregAPA <- calcApa(bedpe = gained,
                         hic = GREG_hic,
                         norm = norm,
                         res = res,
                         buffer = buffer) 

gainedLeukAPA <- calcApa(gained,
                         hic = LEUK_hic,
                         norm = norm,
                         res = res,
                         buffer = buffer) 

lostGregAPA <-calcApa(lost,
                      hic = GREG_hic,
                      norm = norm,
                      res = res,
                      buffer = buffer)

lostLeukAPA <-calcApa(lost,
                      hic = LEUK_hic,
                      norm = norm,
                      res = res,
                      buffer = buffer)

## plot APA plots 
library(plotgardener)
library(RColorBrewer)

# set up plotgardener page
pageCreate(width = 6, height = 5, showGuides = F,
           xgrid = 0, ygrid = 0)

## Define parameters
p <- pgParams(assembly = "hg19",
              x = 0.5,
              y = 0.5,
              width = 1.5,
              height = 1.5,
              space = 0.1,
              palette = colorRampPalette(brewer.pal(9,"YlGnBu")),
              fontcolor = "black")

## APAs normalized to the number of interactions in each category
apas <- list(staticGregAPA/length(static),
             gainedGregAPA/length(gained),
             lostGregAPA/length(lost),
             staticLeukAPA/length(static),
             gainedLeukAPA/length(gained),
             lostLeukAPA/length(lost))

## Specify x and y positions
xpos <- rep(c(p$x, p$x+p$width+p$space, p$x+(p$width+p$space)*2), 2)
ypos <- rep(c(p$y, p$y+p$height+p$space), each = 3)

## Extract center pixel and calculate zrange
cp <- unlist(lapply(apas, \(x) x['0','0']))
zrng <- rep(list(c(0, max(cp[1], cp[4])),
                 c(0, max(cp[2], cp[5])),
                 c(0, max(cp[3], cp[6]))), 2)

## Arrange plots
plots <- 
  lapply(1:6, \(i){
    plotApa(apa = apas[[i]],
            params = p,
            x = xpos[i],
            y = ypos[i],
            zrange = zrng[[i]])
  })

## Add heatmap legends
lapply(4:6, \(i) {
  annoHeatmapLegend(plot = plots[[i]],
                    orientation = 'h',
                    params = p,
                    x = xpos[i],
                    y = ypos[i] + p$height + p$space,
                    height = p$space)
})

## Column labels
plotText(label = c("Static", "PHF23-Loops", "N-IDR_FS/A9-Loops"),
         x = ((((xpos + p$width) - xpos) / 2) + xpos)[1:3],
         y = p$y - p$space,
         just = c('center', 'bottom'),
         fontsize = 8)

## Row labels
plotText(label = c("PHF23", "N-IDR_FS/A9"),
         x = p$x - p$space,
         y = ((((ypos + p$height) - ypos) / 2) + ypos)[c(1,4)],
         just = c('center', 'bottom'),
         rot = 90,
         fontsize = 8)

## Legend label
plotText(label = "Counts per loop",
         x = (((xpos[3] + p$width) - xpos[1])/2) + xpos[1],
         y = ypos[4] + p$height + p$space*3,
         fontsize = 8)




