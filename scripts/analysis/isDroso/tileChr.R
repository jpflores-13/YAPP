## Load libraries
library(plotgardener)
library(GenomicRanges)
library(RColorBrewer)
library(purrr)
library(glue)
library(colorspace)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BiocParallel)

# Set up tiles ------------------------------------------------------------

tile <- 
  tileGenome(seqlengths = seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene),
                   tilewidth = 25e5,
                   cut.last.tile.in.chrom = T) |> 
  keepStandardChromosomes(species = "Homo_sapiens", pruning.mode = "coarse")

tile <- tile[seqnames(tile) == "chr1"]

# Set-up params -----------------------------------------------------------
pdf(file = "plots/tile_chr1.pdf",
    width = 5.5,
    height = 7.5)

for (i in seq_along(tile)){
  
  p <- pgParams(assembly = "hg38",
                resolution = 5e3,
                chrom = gsub('chr', '',
                             as.character(seqnames(tile)[i])),
                chromstart = start(tile)[i],
                chromend = end(tile)[i],
                zrange = c(0,100),
                norm = "KR",
                x = 0.25,
                width = 5,
                length = 5,
                height = 2,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin visualization -----------------------------------------------------
  
  ## Make page
  pageCreate(width = 5.5, height = 7.5,
             xgrid = 0, ygrid = 0, showGuides = F)

  
  ## Plot Hi-C
  omega <- plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic",
                            params = p,
                            y = 0.5)
  
  cont <- plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                           params = p,
                           y = 2.5)
  
  sorb <- plotHicRectangle(data = "data/raw/hic/hg38/220628_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic",
                           params = p,
                           y = 4.5)
  
  ## Plot Gene Track
  plotGenes(params = p,
            chrom = paste0('chr', p$chrom),
            height = 0.75,
            y = 6.5)
  
  
  ## Plot Gene Region
  plotGenomeLabel(params = p,
                  y = 7.25)
  
  ## Plot Text
  plotText(label = "omega",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"))
  
  plotText(label = "control",
           x = 0.25,
           y = 2.5,
           just = c("top", "left"))
  
  plotText(label = "sorbitol",
           x = 0.25,
           y = 4.5,
           just = c("top", "left"))
}

dev.off()