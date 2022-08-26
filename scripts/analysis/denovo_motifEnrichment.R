## De novo (unbiased) motif enrichment at gained and lost loop anchors

# Creating focal and background ATAC peaks --------------------------------

## library set-up
library(InteractionSet)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(nullranges)
library(RColorBrewer)
library(glue)
library(hictoolsr)
library(memes)
library(plyranges)
library(tidyverse)
library(glue)

# Data Import -------------------------------------------------------------
## loops
loops <- readRDS("data/processed/hic/YAPP_hic_loopCounts.rds")

## add a loop_count column
loops$loop_name <- glue("loop_{1:length(loops)}")

## ATAC peaks
peaks <- read.table("data/raw/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv", header = TRUE) 

# Background Peak Set-up --------------------------------------------------
## all ATAC peaks at ANY loop anchor
background_peaks <- GRanges(seqnames = Rle(peaks$chr), 
                            ranges = IRanges(start = peaks$start, end = peaks$stop), 
                            peak = rownames(peaks))

# Foreground Peak set-up --------------------------------------------------
## import differential loops
diffLoops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

## all ATAC peaks at gained loop anchors
gained_loops <- diffLoops |> 
  subset(padj < 0.1 & log2FoldChange > 0)

## all ATAC peaks at lost loop anchors
lost_loops <- diffLoops |> 
  subset(padj < 0.1 & log2FoldChange < 0)

focal_peaks_gained  <- subsetByOverlaps(background_peaks, gained_loops)
focal_peaks_lost <- subsetByOverlaps(background_peaks, lost_loops)

# Prep for MEME ----------------------------------------------------------
human.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

## for MEME, need to make sure you're using 1-base ranges
background_peaks[seqnames(background_peaks) == "chr9_KI270718v1_random"] <- 
  shift(background_peaks[seqnames(background_peaks) == "chr9_KI270718v1_random"], 1)  

background_peaks[seqnames(background_peaks) == "chrUn_GL000195v1"] <- 
  shift(background_peaks[seqnames(background_peaks) == "chrUn_GL000195v1"], 1)  

background_peaks[seqnames(background_peaks) == "chrUn_KI270333v1"] <- 
  shift(background_peaks[seqnames(background_peaks) == "chrUn_KI270333v1"], 1)  

background_peaks[seqnames(background_peaks) == "chrUn_KI270336v1"] <- 
  shift(resize(background_peaks[seqnames(background_peaks) == "chrUn_KI270336v1"],
               width = length(human.genome$chrUn_KI270336v1)), 1)

background_peaks[seqnames(background_peaks) == "chrUn_KI270337v1"] <- 
  shift(background_peaks[seqnames(background_peaks) == "chrUn_KI270337v1"], 1) 

background_peaks[seqnames(background_peaks) == "chrUn_KI270442v1"] <- 
  shift(background_peaks[seqnames(background_peaks) == "chrUn_KI270442v1"], 1)

background_peaks[seqnames(background_peaks) == "chrUn_KI270442v1"][41] <- 
  resize(background_peaks[seqnames(background_peaks) == "chrUn_KI270442v1"][41], width = 1041)

## get sequences of background peaks
sequence_background <- background_peaks |>
  get_sequence(human.genome)

## get sequences of focal peaks
sequence_focal_gained <- focal_peaks_gained |> 
  get_sequence(human.genome) 

sequence_focal_lost <- focal_peaks_lost |> 
  get_sequence(human.genome) 

options(meme_db = "data/raw/atac/meme_files/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme")

de_novo_gained <- runStreme(input = sequence_focal_gained, control = sequence_background,
                            outdir = "tables/atac/de_novo/gained")

View(de_novo_gained)

de_novo_lost <- runStreme(input = sequence_focal_lost, control = sequence_background,
                          outdir = "tables/atac/de_novo/lost")

View(de_novo_lost)

Visualization -----------------------------------------------------------
  de_novo_gained_top <- de_novo_gained |> 
  mutate(log10pval = (-log10(pvalue))) |>
  mutate(motif_id = str_remove(motif_id, "_.*")) |> 
  arrange(rank) |> 
  slice_head(n = 50)

de_novo_lost_top <- de_novo_lost |> 
  mutate(log10pval = (-log10(pvalue))) |> 
  mutate(motif_id = str_remove(motif_id, "_.*")) |> 
  arrange(rank) |> 
  slice_head(n = 50)

top_de_novo_gained <- de_novo_gained_top |> 
  ggplot(aes(x = reorder(motif_id, -log10pval), y = log10pval)) +
  geom_col(fill = "steelblue") +
  labs(title = "Motif Enrichment at Gained Loop Anchors", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))

top_de_novo_lost <- de_novo_lost_top |> 
  ggplot(aes(x = reorder(motif_id, -log10pval), y = log10pval)) +
  geom_col(fill = "steelblue") +
  labs(title = "Motif Enrichment at Lost Loop Anchors", x = "motif_id") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))

library(plotgardener)

##make pdf
pdf(file = "plots/denovo_motifEnrichment_meme.pdf",
    width = 10,
    height = 7)

# Begin Visualization -----------------------------------------------------
## Make page
pageCreate(width = 10, height = 7,
           xgrid = 0, ygrid = 0, showGuides = F)

plotGG(top_de_novo_gained, 
       x = 0.5,
       y = 0.5,
       width = 9,
       height = 3)

plotGG(top_de_novo_lost,
       x = 0.5,
       y = 3.5,
       width = 9, 
       height = 3)

dev.off()








