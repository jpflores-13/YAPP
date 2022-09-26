# Motif Enrichment Analysis using monaLisa

## library set-up
library(SummarizedExperiment)
library(JASPAR2020)
library(TFBSTools)
library(monaLisa)
library(ComplexHeatmap)
library(circlize)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(nullranges)
library(RColorBrewer)
library(hictoolsr)
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

# Prep for  ----------------------------------------------------------
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

# Prep for HOMER ----------------------------------------------------------
txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                            release=NA,
                            circ_seqs=NULL,
                            server="ensembldb.ensembl.org",
                            username="anonymous", password=NULL, port=0L,
                            tx_attrib=NULL)

genes <- genes(txdb)
genes <- data.frame(genes)
genes$seqnames <- paste0("chr", genes$seqnames)
genes <- GRanges(seqnames = Rle(genes$seqnames), 
                 ranges = IRanges(start = genes$start, end = genes$end), 
                 strand = genes$strand,
                 ensembl_gene_id = genes$gene_id)


# HOMER/ monaLisa ---------------------------------------------------------

# bin regions
# - peak_change is a numerical vector
# - peak_change needs to be created by the user to run this code
peak_bins <- bin(x = peak_change, binmode = "equalN", nElement = 400)

# calculate motif enrichments
# - peak_seqs is a DNAStringSet, pwms is a PWMatrixList
# - peak_seqs and pwms need to be created by the user to run this code
se <- calcBinnedMotifEnrR(seqs = peak_seqs,
                          bins = peak_bins,
                          pwmL = pwms)


