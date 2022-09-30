## How many gained ATAC regions overlap with gene promoters?
## How many static/gained/loop anchors overlap with gene promoters?

## Load libraries -----------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(InteractionSet)
library(nullranges)
library(RColorBrewer)
library(ggplot2)

# data import -------------------------------------------------------------

ATAC <- readRDS("data/processed/atac/YAPP_hic_diff_ATACcounts.rds")
gainedATAC <- ATAC |> 
  subset(log2FoldChange > 0 & padj < 0.1)

loops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

## separate loops into gained, lost, static
gained <- loops |>
  subset(padj < 0.1 & log2FoldChange > 0)

lost <- loops |>
  subset(padj < 0.1 & log2FoldChange < 0)

static <- loops |>
  subset(padj > 0.1)


# gene promoter set-up ----------------------------------------------------


#Genes
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
promoters <- promoters(genes)

promoter_gainedATAC_overlap <- 
  subsetByOverlaps(promoters, gainedATAC) |> 
  data.frame()


# gained atac regions + promoters -----------------------------------------

library(tidyverse)  
promoter_gainedATAC_overlap <- promoter_gainedATAC_overlap |> 
  dplyr::select(ensembl_gene_id)

write.csv(promoter_gainedATAC_overlap, file = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/tables/atac/de_novo/gained/streme/promoter_gainedATAC_overlap.csv")

# loop anchor + promoter overlaps -----------------------------------------

gainedLoop_promoter_overlaps <- subsetByOverlaps(promoters, gained)
write.csv(gainedLoop_promoter_overlaps, file = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/tables/hic/gainedLoop_promoter_overlaps.csv")

staticLoop_promoter_overlaps <- subsetByOverlaps(promoters, static)
write.csv(staticLoop_promoter_overlaps, file = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/tables/hic/staticLoop_promoter_overlaps.csv")

lostLoop_promoter_overlaps <- subsetByOverlaps(promoters, lost)
write.csv(lostLoop_promoter_overlaps, file = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/tables/hic/lostLoop_promoter_overlaps.csv")


