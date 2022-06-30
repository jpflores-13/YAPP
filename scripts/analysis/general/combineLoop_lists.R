## Combine Loops Lists

## install github
# remotes::install_github("EricSDavis/hictoolsr")

## load package
library(hictoolsr)
library(data.table)
library(tidyverse)
library(InteractionSet)
library(dbscan)

# for merging loops -----------------------------
loop_files <- c("data/hic/isDroso/sorb/loops/5kbLoops.txt",
                "data/hic/isDroso/cont/loops/5kbLoops.txt",
                "data/hic/noDroso/cont/loops/5kbLoops.txt",
                "data/hic/noDroso/sorb/loops/5kbLoops.txt",
                "data/hic/isDroso/cont/loops/5kbLoops-omega.txt")

# for extracting counts -----------------------------
hic_files <- c("data/hic/isDroso/cont/YAPP_HEK_control_1_inter_30.hic",
               "data/hic/isDroso/cont/YAPP_HEK_control_2_inter_30.hic",
               "data/hic/isDroso/cont/YAPP_HEK_control_3_inter_30.hic",
               "data/hic/isDroso/sorb/YAPP_HEK_sorbitol_4_inter_30.hic",
               "data/hic/isDroso/sorb/YAPP_HEK_sorbitol_5_inter_30.hic",
               "data/hic/isDroso/sorb/YAPP_HEK_sorbitol_6_inter_30.hic")

# merge / extract -----------------------------
loopCounts <- mergeBedpe(bedpeFiles = loop_files,
                         res = 10000,
                         selectCol = 12,
                         dist_method = "manhattan",
                         minPts = 2) |> 
  as_ginteractions() |> 
  binBedpe(res = 10e3, 
           a1Pos = "center", 
           a2Pos = "center") |> 
  
  extractCounts(hic = hic_files,
                chroms = c(1:22, "X", "Y"),
                res = 10e3,
                norm = "NONE",
                matrix = "observed")

# save data -----------------------------
save(loopCounts, file = "data/output/isDroso/loopCounts.rda")
