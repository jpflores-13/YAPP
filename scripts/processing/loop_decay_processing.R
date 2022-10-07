library(strawr)
library(dbscan)
library(InteractionSet)
library(raster)
library(hictoolsr)
library(glue)
library(hictoolsr)
library(mariner)
library(plyranges)
library(nullranges)
library(tidyverse)

## load function from utils folder
source("scripts/utils/mh_index.R")

## loop upload
all_loops <- readRDS("data/processed/hic/YAPP_hic_diff_loopCounts.rds")

gained_loops <- all_loops |> 
  subset(padj < 0.1 & log2FoldChange > 0)

lost_loops <- all_loops |> 
  subset(padj < 0.1 & log2FoldChange < 0)

## create metadata columns for contact frequency & size 
mcols(all_loops)$loop_size <- pairdist(all_loops)
mcols(all_loops)$loop_type <- case_when(
  mcols(all_loops)$padj < 0.1 & mcols(all_loops)$log2FoldChange > 0 ~ "gained",
  mcols(all_loops)$padj < 0.1 & mcols(all_loops)$log2FoldChange < 0 ~ "lost",
  mcols(all_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

## contact frequency
mcols(all_loops)$sorb_contacts <- mcols(all_loops)$YAPP_HEK_sorbitol_4_2_inter_30.hic +
  mcols(all_loops)$YAPP_HEK_sorbitol_5_2_inter_30.hic + mcols(all_loops)$YAPP_HEK_sorbitol_6_2_inter_30.hic

## use matchRanges to select a null set of control sample loops that is matched for size & contact frequency
nullSet <- matchRanges(focal = all_loops[mcols(all_loops)$loop_type == "gained"],
                       pool = all_loops[!mcols(all_loops)$loop_type == "gained"],
                       covar = ~ loop_size + sorb_contacts, 
                       method = 'stratified',
                       replace = FALSE)

## bring in merged sorbitol .hic file 
sorb_hic <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic"
cont_hic <- "data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic"

# loop decay for gained loops within sorbitol .hic file ---------------------------------------

# calculate APA matrices for sorb .hic file
normalized_gained=data.frame()
for(j in 1:length(gained_loops)){
  l <- calcApa(bedpe = gained_loops [j], hic = sorb_hic, norm = "NONE", res = 10e3, buffer = 10, 
               filter = FALSE) 
  for(i in 0:10){
    normalized_gained[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}

head(normalized_gained)
colnames(normalized_gained) <- c("0","1","2","3","4","5","6","7","8","9","10")

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero_gained=data.frame()
for(i in 1:length(gained_loops)){
  for(j in 1:11){
    a <- normalized_gained[i,1]-normalized_gained[i,11]
    b <- normalized_gained[i,j]-normalized_gained[i,11]
    zero_gained[i,j] <- round((b/a),2)
  }
  
}
head(zero_gained)
colnames(zero_gained) <- c("0","1","2","3","4","5","6","7","8","9","10")


Group_gained <- factor(0:10)
normalized1_gained <- na.omit(zero_gained)
saveRDS(normalized1_gained, "data/processed/hic/gained_normalized_mh_index.rds")

# loop decay for lost loops within sorbitol .hic file ---------------------------------------

# calculate APA matrices for sorb .hic file
normalized_lost=data.frame()
for(j in 1:length(lost_loops)){
  l <- calcApa(bedpe = lost_loops [j], hic = sorb_hic, norm = "NONE", res = 10000, buffer = 10, 
               filter = FALSE) 
  for(i in 0:10){
    normalized_lost[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}

head(normalized_lost)
colnames(normalized_lost) <- c("0","1","2","3","4","5","6","7","8","9","10")

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero_lost=data.frame()
for(i in 1:length(lost_loops)){
  for(j in 1:11){
    a <- normalized_lost[i,1]-normalized_lost[i,11]
    b <- normalized_lost[i,j]-normalized_lost[i,11]
    zero_lost[i,j] <- round((b/a),2)
  }
}
head(zero_lost)
colnames(zero_lost) <- c("0","1","2","3","4","5","6","7","8","9","10")

Group_lost <- factor(0:10)
normalized1_lost <- na.omit(zero_lost)
saveRDS(normalized1_lost, "data/processed/hic/lost_normalized_mh_index.rds")

# loop decay for all loops within control .hic file ---------------------------------------

# calculate APA matrices for sorb .hic file
normalized_control=data.frame()
for(j in 1:length(all_loops)){
  l <- calcApa(bedpe = all_loops [j], hic = cont_hic, norm = "NONE", res = 10000, buffer = 10, 
               filter = FALSE) 
  for(i in 0:10){
    normalized_control[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}

head(normalized_control)
colnames(normalized_control) <- c("0","1","2","3","4","5","6","7","8","9","10")

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero_control=data.frame()
for(i in 1:length(all_loops)){
  for(j in 1:11){
    a <- normalized_control[i,1]-normalized_control[i,11]
    b <- normalized_control[i,j]-normalized_control[i,11]
    zero_control[i,j] <- round((b/a),2)
  }
}
head(zero_control)
colnames(zero_control) <- c("0","1","2","3","4","5","6","7","8","9","10")


Group_control<- factor(0:10)
normalized1_control <- na.omit(zero_control)
saveRDS(normalized1_control, "data/processed/hic/cont_normalized_mh_index.rds")

# loop decay for nullSet loops within control .hic file ---------------------------------------
## MatchedGInteractions to GInteractions
nullSet <- nullSet |>
  data.frame() |>
  mariner::as_ginteractions()

# calculate APA matrices for sorb .hic file
normalized_nullSet=data.frame()
for(j in 1:length(nullSet)){
  l <- calcApa(bedpe = nullSet [j], hic = cont_hic, norm = "NONE", res = 10000, buffer = 10, 
               filter = FALSE) 
  for(i in 0:10){
    normalized_nullSet[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}

head(normalized_nullSet)
colnames(normalized_nullSet) <- c("0","1","2","3","4","5","6","7","8","9","10")

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero_nullSet=data.frame()
for(i in 1:length(nullSet)){
  for(j in 1:11){
    a <- normalized_nullSet[i,1]-normalized_nullSet[i,11]
    b <- normalized_nullSet[i,j]-normalized_nullSet[i,11]
    zero_nullSet[i,j] <- round((b/a),2)
  }
}
head(zero_nullSet)
colnames(zero_nullSet) <- c("0","1","2","3","4","5","6","7","8","9","10")


Group_nullSet<- factor(0:10)
normalized1_nullSet <- na.omit(zero_nullSet)
saveRDS(normalized1_nullSet, "data/processed/hic/nullSet_cont_normalized_mh_index.rds")
