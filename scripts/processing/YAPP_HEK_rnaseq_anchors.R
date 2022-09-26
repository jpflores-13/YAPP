## RNA-seq Differential Gene Expression Analysis for YAPP HEK Hi-C static, lost, and gained loops

## Load required libraries
library(data.table)
library(plyranges)
library(tximeta)
library(readr)
library(glue)
library(DESeq2)
library(InteractionSet)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Read in sample sheet as colData
colData <-
  fread("data/raw/rna/output/YAPP_HEK_WT_1_RNApipeSamplesheet.txt") |> 
  as.data.frame()

## Use 0hr and 6hr time points 
colData <- colData |> 
  dplyr::filter(Time %in% c("0h", "1h"))
  
## Edit quant paths 
colData$files <- paste0("data/raw/rna/output/quant/", colData$sn, "/quant.sf")
setnames(colData, "sn", "names")

## Check that quant paths are correct
file.exists(colData$files)

## Import data with tximeta & summarize to gene
se <- tximeta(colData)
gse <- summarizeToGene(se)

## Convert to factors
colData(gse)[] <- lapply(colData(gse), as.factor)

## Build DESeq object
dds <- DESeqDataSet(gse, design = ~Bio_Rep + Treatment)

## Filter out lowly expressed genes (at least 10 counts in at least 4 samples)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]

## Fit model
dds <- DESeq(dds)
resultsNames(dds)

## Get DE gene results as GRanges
de_genes <- results(dds,
                    name = "Treatment_sorb_vs_cont",
                    format = "GRanges") %>%
  names_to_column("gene_id")

## Add gene symbols to de_genes and
## remove non-standard chromosomes
de_genes <- 
  inner_join(x = as.data.table(de_genes),
             y = as.data.table(rowData(gse)) |>
               dplyr::select(c("gene_id", "symbol", "tx_ids")),
             by = "gene_id") |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
  keepStandardChromosomes(pruning.mode = "coarse")

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene |>
  keepStandardChromosomes()

seqlevels(de_genes) <- seqlevels(txdb)
seqinfo(de_genes) <- seqinfo(txdb)

## Save results to a file
saveRDS(de_genes, file = "data/processed/rna/YAPP_HEK_rnaseq_anchors.rds")
