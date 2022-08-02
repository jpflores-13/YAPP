
## use DESeq2 to call differential loops 

library(DESeq2)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
library(pheatmap)
library(apeglm)

# Load Data--------------
load("data/processed/hic/YAPP_hic_loopCounts.rda")

# Add a loop_count column --------------
loopCounts$loop_name <- glue("loop_{1:length(loopCounts)}")

# Create a matrix --------------
m <- mcols(loopCounts)[, grep("*inter.*", colnames(mcols(loopCounts)))] %>% 
  as.matrix()

# Construct colData/metadata ----------------------------------------------

## String split the colnames 

colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")), stringsAsFactors = T)
rownames(colData) <- colnames(m)
colnames(colData) <- c("Project", "Cell_Type", "Treatment", "Replicate")
colData <- colData[, c(1:4)]

## Correct replicate number
colData <- colData |> 
  mutate(Replicate = str_replace(Replicate, "4", "1")) |> 
  mutate(Replicate = str_replace(Replicate, "5", "2")) |> 
  mutate(Replicate = str_replace(Replicate, "6", "3"))

## Make sure sample names match
all(colnames(m) == rownames(colData))

# Run DESeq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = m,
                       colData = colData,
                       design = ~Replicate + Treatment)

## Hypothesis testing with Wald with `betaPrior = F`. No LRT because no time points in this analysis.
dds <- DESeq(dds)

## Plot dispersion estimates
plotDispEsts(dds)

# QC via data visualization before moving on ------------------------------

## plot PCA
pdf(file = "plots/YAPP_diffLoops_sorb_PCA_hic.pdf")
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
pdf(file = "plots/YAPP_diffLoops_sorb_hcluster_hic.pdf")
pheatmap(rld_cor)
dev.off()

## outputs
res <- results(dds)
resultsNames(dds)
summary(res)
res <- lfcShrink(dds, coef="Treatment_sorbitol_vs_control", type= "apeglm")

summary(res)
pdf(file = "plots/YAPP_diffLoops_sorb_MA_hic.pdf")
plotMA(res, ylim=c(-4,4), main = "Differential Loops - Sorbitol-Treated (top) vs. Untreated (bottom)",
       ylab = "LFC",
       xlab = "mean of norm. counts")
dev.off()

# Concatenate loopCounts and DESeqResults ---------------------------------
mcols(loopCounts) <- cbind(mcols(loopCounts), res)
diff_loopCounts <- loopCounts

## save as .rda
save(diff_loopCounts, file = "data/processed/hic/YAPP_hic_diff_loopCounts.rda")
sessionInfo()
