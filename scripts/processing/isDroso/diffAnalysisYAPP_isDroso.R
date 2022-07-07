
## use DESeq2 to call differential loops 

library(DESeq2)
library(dplyr)
library(glue)
library(stringr)
library(purrr)

## construct your matrix
#rows = loops
#cols = samples

# Load Data--------------
load("data/processed/hic/isDroso/loopCounts.rda")

# Add a loop_count column --------------
loopCounts$loop_name <- glue("loop_", "{1:length(loopCounts)}")

# Create a matrix --------------
m <- mcols(loopCounts)[, grep("*inter.*", colnames(mcols(loopCounts)))] %>% 
  as.matrix()

# Change colnames  --------------
rep <- rep(1:3, times = 2)
cond <- rep(c("control", "sorbitol"), each = 3)
colnames(m) <- c(glue("YAPP_HEK_{cond}_{rep}_inter_30.hic"))

# Assign rownames  --------------
rownames(m) <- glue("{loopCounts$loop_name}")

#construct your colData  --------------

# String split the colnames --------------

colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")), stringsAsFactors = T)

rownames(colData) <- colnames(m)

## colData should be a dataframe
## colnames(m) should be your rownames 

colData <- colData[, c(1:4)]

colnames(colData) <- c("Project", "Cell_Type", "Treatment", "Replicate")

## DESeq2 analysis

dds <- DESeqDataSetFromMatrix(countData = m,
                              colData = colData,
                              design = ~ Replicate + Treatment)

## Trying different statistical tests (LRT vs. Wald)

# dds <- DESeq(dds, test = "LRT", reduced = ~ Replicate)
dds <- DESeq(dds, test = "Wald", betaPrior = T)
# dds <- DESeq(dds, test = "Wald")
res <- results(dds)

resultsNames(dds)
# res <- lfcShrink(dds, coef="Treatment_sorbitol_vs_control", type= "apeglm")

summary(res)
pdf(file = "plots/YAPP_diffLoops_sorb_MA_isDroso_betaPriorT.pdf")
plotMA(res, ylim=c(-4,4), main = "Differential Loops - Untreated vs. Sorbitol-Treated",
       ylab = "LFC",
       xlab = "mean of norm. counts")
dev.off()

plotPCA(vst(dds), intgroup = "Treatment") + ggplot2::theme(aspect.ratio = 1)
plotPCA(vst(dds), intgroup = "Treatment", returnData = TRUE)

#concatenate loopCounts and DESeqResults
all(rownames(res) == loopCounts$loop_name)
mcols(loopCounts) <- cbind(mcols(loopCounts), res)
diff_loopCounts <- loopCounts

save(diff_loopCounts, file = "data/output/hic/isDroso/diff_loopCounts.rda")

