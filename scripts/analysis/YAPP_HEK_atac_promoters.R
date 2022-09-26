## Load libraries -----------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(InteractionSet)
library(nullranges)
library(RColorBrewer)
library(ggplot2)


## Read in data ----------------------------------------------------------------------

# Read in text files containing either all or differential elements 
#Loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> 
  as.data.frame()
loops <- loops[,c(1:3,6:8)]
rownames(loops) <- paste0("loop", 1:nrow(loops))

#ATAC
allATAC <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacLFC.rds")


## Convert to GRanges: 
#Loops
loopAnc <- GInteractions(anchor1 = GRanges(seqnames = Rle(loops$seqnames1), 
                                           ranges = IRanges(start = loops$start1, end = loops$end1)), 
                         anchor2 = GRanges(seqnames = Rle(loops$seqnames2), 
                                           ranges = IRanges(start = loops$start2, end = loops$end2)), 
                         name = rownames(loops))
#ATAC 
allATACGR <- GRanges(seqnames = Rle(allATAC$chr), 
                     ranges = IRanges(start = allATAC$start, end = allATAC$stop), 
                     peak = rownames(allATAC))

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
staticPromoters <- promoters[!promoters$ensembl_gene_id %in% diffGenes,]

#Get the coordinate information for each of the genes in each of the clusters: 
for (i in 1:length(clusters)){
  c <- clusters[[i]]
  clusters[[i]] <- promoters[promoters$ensembl_gene_id %in% c$V1,] 
}