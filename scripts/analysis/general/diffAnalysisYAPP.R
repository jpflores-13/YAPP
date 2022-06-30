
## use DESeq2 to call differential loops 

library(DESeq2)
library(tidyverse)

## construct your matrix
#rows = loops
#cols = samples


# Load Data--------------
load("data/output/loopCounts.rda")

# Add a loop_count column --------------
loopCounts$loop_name <- paste0("loop_", 1:length(loopCounts))

# Create a matrix --------------
m <- mcols(loopCounts)[, grep("*inter.*", colnames(mcols(loopCounts)))] %>% 
  as.matrix()

# Change colnames  --------------
colnames(m) <- c("YAPP_HEK_control_1_inter_30.hic",
                 "YAPP_HEK_control_2_inter_30.hic",
                 "YAPP_HEK_control_3_inter_30.hic",
                 "YAPP_HEK_sorbitol_4_inter_30.hic",
                 "YAPP_HEK_sorbitol_5_inter_30.hic",
                 "YAPP_HEK_sorbitol_6_inter_30.hic")

# Assign rownames  --------------
rownames(m) <- paste0(loopCounts$loop_name)

#construct your colData  --------------

# String split the colnames --------------

colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")))

rownames(colData) <- colnames(m)

## colData should be a dataframe
## colnames(m) should be your rownames 

colData <- colData[, c(1, 2, 3, 4)]

colnames(colData) <- c("Project", "Cell_Type", "Treatment", "Replicate")

# Specify Accurate Replicate Numbers  --------------
colData <- colData |> 
  mutate(Replicate = str_replace(Replicate, "4", "1")) |> 
  mutate(Replicate = str_replace(Replicate, "5", "2")) |> 
  mutate(Replicate = str_replace(Replicate, "6", "3"))
  

## DESeq2 analysis

dds <- DESeqDataSetFromMatrix(countData = m,
                              colData = colData,
                              design = ~ Replicate + Treatment)

dds <- DESeq(dds)
res <- results(dds)

resultsNames(dds)
res <- lfcShrink(dds, coef="Treatment_sorbitol_vs_control", type= "apeglm")
# From Andrea ---------------------
# p = 0.05
# L = 0
# sig_LRT <- res[which(res$padj < p & abs(res$log2FoldChange) >= L),]

summary(res)

plotMA(res, ylim=c(-4,4), main = "Differential Loops - Untreated vs. Sorbitol-Treated",
       ylab = "LFC",
       xlab = "mean of norm. counts")
# plotPCA(vst(dds), intgroup = "Treatment") + ggplot2::theme(aspect.ratio = 1)
plotPCA(vst(dds), intgroup = "Treatment", returnData = TRUE)

#concatenate loopCounts and DESeqResults
all(rownames(res) == loopCounts$loop_name)
mcols(loopCounts) <- cbind(mcols(loopCounts), res)
diff_loopCounts <- loopCounts

save(diff_loopCounts, file = "data/output/diff_loopCounts.rda")

