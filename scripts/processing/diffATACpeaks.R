## Perform differential loops analysis with DESeq2 with sizeFactors = 1

library(DESeq2)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
library(pheatmap)
library(apeglm)

# Load ATAC peak table ---------------------------------------------

peaks <- read.table("data/raw/atac/output/peaks/YAPP_HEK_1_peakCounts.tsv",
                    header = T)

## add a peak_count column
peaks$peak_count <- paste0("peak_", 1:nrow(peaks))

## convert to GRanges object
peaks <- GRanges(seqnames = Rle(peaks$chr), 
        ranges = IRanges(start = peaks$start, end = peaks$stop),
        YAPP_HEK_cont_0h_1_1 = peaks[,4],
        YAPP_HEK_cont_0h_2_1 = peaks [,5],
        YAPP_HEK_cont_0h_3_1 = peaks[,6],
        YAPP_HEK_sorb_1h_2_1 = peaks[,7],
        YAPP_HEK_sorb_1h_3_1 = peaks[,8],
        YAPP_HEK_sorb_1h_4_1 = peaks[,9],
        peak_count = peaks[,10])

# Create matrices for countData -------------------------------------------

m <- mcols(peaks)[, grep("YAPP.*", colnames(mcols(peaks)))] %>% 
  as.matrix()

# Construct colData/metadata ----------------------------------------------

## String split the colnames 

colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")), stringsAsFactors = T)
rownames(colData) <- colnames(m)
colnames(colData) <- c("Project", "Cell_Type", "Treatment", "Time", "Replicate")
colData <- colData[, c(1:3,5)]

## Correct replicate number
colData <- colData |> 
  mutate(Replicate = rep(c(1:3), 2))

colData$Replicate <- as.factor(colData$Replicate)

## Make sure sample names match
all(colnames(m) == rownames(colData))

# Run DESeq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = m,
                              colData = colData,
                              design = ~ Replicate + Treatment)

## disable DESeq's default normalization 
# sizeFactors(dds) <- rep(1, ncol(dds))

## Hypothesis testing with Wald with `betaPrior = F`
dds <- DESeq(dds)

## Plot dispersion estimates
plotDispEsts(dds)

# QC via data visualization before moving on ------------------------------

## plot PCA
pdf(file = "plots/YAPP_diffATAC_sorb_PCA.pdf")
plotPCA(vst(dds), intgroup = "Treatment") + ggplot2::theme(aspect.ratio = 1)
dev.off()
plotPCA(vst(dds), intgroup = "Treatment", returnData = TRUE)

## transform counts for hierarchical clustering
rld <- rlog(dds, blind=TRUE)

## Extract the rlog matrix from the object
rld_mat <- assay(rld)

## Compute pairwise correlation values
rld_cor <- cor(rld_mat)

## Plot heatmap
pdf(file = "plots/YAPP_diffATAC_sorb_hcluster.pdf")
pheatmap(rld_cor)
dev.off()

## outputs
res <- results(dds)
resultsNames(dds)
summary(res)
res <- lfcShrink(dds, coef="Treatment_sorb_vs_cont", type= "apeglm")

summary(res)
pdf(file = "plots/YAPP_diffATAC_sorb_MA.pdf")
plotMA(res, ylim=c(-4,4), main = "Differential ATAC Analysis",
       ylab = "LFC",
       xlab = "mean of norm. counts")
dev.off()

# Concatenate loopCounts and DESeqResults ---------------------------------
mcols(peaks) <- cbind(mcols(peaks), res)
diff_peakCounts <- peaks

## save as .rds
saveRDS(peaks, file = "data/processed/atac/YAPP_hic_diff_ATACcounts.rds")
sessionInfo()
