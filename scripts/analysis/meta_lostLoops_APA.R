## Lost Loop APA plots for HEK Hi-C, HEK Micro-C, and T47D Hi-C data

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

## Filter out lost control loops
lostLoops <- diff_loopCounts |> 
  subset(padj < 0.1 & log2FoldChange < 0)

## Define resolution, buffer, and normalization (pixels from center)
res <- 10e3
buffer <- 10
norm <- "SCALE"

## Filter out short loop interactions
filteredLoops <- lostLoops |> 
  filterBedpe(res = res, 
              buffer = buffer) |> 
  'names<-'(value = names(lostLoops))

## Summary loops
filteredLoops |> 
  summary()

## Hi-C file paths 
cont_hic_YAPP <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic"
sorb_hic_YAPP <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic"
cont_hic_HYPE <- "external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic"
nacl_hic_HYPE <- "external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic"

## Calculate APA matrices for loops from control HEK Hi-C data
cont_APA_mat_hic_YAPP <- filteredLoops |> 
  calcApa(hic = cont_hic_YAPP, norm = norm, buffer = buffer)

## Calculate APA matrices for loops from sorbitol HEK Hi-C data
sorb_APA_mat_hic_YAPP <- filteredLoops |> 
  calcApa(hic = sorb_hic_YAPP, norm = norm, buffer = buffer)

## Calculate APA matrices for loops from control T47D Hi-C data
cont_APA_mat_hic_HYPE <- filteredLoops |> 
  calcApa(hic = cont_hic_HYPE, norm = norm, buffer = buffer)

## Calculate APA matrices for loops from NaCl T47D Hi-C data
nacl_APA_mat_hic_HYPE <-filteredLoops |> 
  calcApa(hic = nacl_hic_HYPE, norm = norm, buffer = buffer)

## Get the number of loops for each condition
nLoops <- filteredLoops |> 
  length()

## Divide each matrix by nLoops
cont_APA_mat_hic_YAPP <- (cont_APA_mat_hic_YAPP/nLoops)
sorb_APA_mat_hic_YAPP <- (sorb_APA_mat_hic_YAPP/nLoops)
cont_APA_mat_hic_HYPE <- (cont_APA_mat_hic_HYPE/nLoops)
nacl_APA_mat_hic_HYPE <- (nacl_APA_mat_hic_HYPE/nLoops)


## combine APA matrices to pull out the max value for zrange max
mats_combined <- c(cont_APA_mat_hic_YAPP,
                   sorb_APA_mat_hic_YAPP,
                   cont_APA_mat_hic_HYPE,
                   nacl_APA_mat_hic_HYPE)


# Create plotgardener page ------------------------------------------------
pdf(file = "plots/meta_lostLoops_APA.pdf",
    height = 3,
    width = 4.25)
a <- dev.cur()

png(file = "vignettes/assets/meta_lostLoops_APA.png",
    height = 3,
    width = 4.25)
dev.control("enable")

## Initiate plotgardener page
pageCreate(width = 3, height = 3, showGuides = F)

## Define shared parameters
p <- pgParams(x = 0.5,
              y = 0.5,
              width = 1,
              height = 1,
              space = 0.075,
              zrange = c(0, max(mats_combined)),
              palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))

## plot untreated HEK Hi-C APA matrix
cont_YAPP_APA <- plotApa(params = p, apa = cont_APA_mat_hic_YAPP)

## plot untreated T47D Hi-C APA matrix
cont_HYPE_APA <- plotApa(params = p, apa = cont_APA_mat_hic_HYPE, x = 1.6)

## plot sorbitol-treated HEK Hi-C APA matrix
sorb_apa <- plotApa(params = p, apa = sorb_APA_mat_hic_YAPP, y = 1.6)

## plot NaCl-treated T47D Hi-C APA matrix
nacl_apa <- plotApa(params = p, apa = nacl_APA_mat_hic_HYPE, y = 1.6, x = 1.6)

# Add legend
annoHeatmapLegend(plot = cont_YAPP_APA,
                  orientation = "h",
                  x = 0.5,
                  y = 2.75,
                  width = 1,
                  height = 0.1,
                  fontcolor = 'black')

# Add legend
annoHeatmapLegend(plot = cont_HYPE_APA,
                  orientation = "h",
                  x = 1.65,
                  y = 2.75,
                  width = 1,
                  height = 0.1,
                  fontcolor = 'black')

## Add text labels
plotText(label = "HEK Hi-C",
         x = 1,
         y = 0.4,
         just = c('center', 'bottom'))

plotText(label = "T47D Hi-C",
         x = 2.1,
         y = 0.4,
         just = c('center', 'bottom'))

plotText(label = "Cont",
         x = 0.4,
         y = 1,
         rot = 90,
         just = c('center', 'bottom'))  

plotText(label = "Treated",
         x = 0.4,
         y = 2.15,
         rot = 90,
         just = c('center', 'bottom'))  
dev.copy(which = a)
dev.off()
dev.off()