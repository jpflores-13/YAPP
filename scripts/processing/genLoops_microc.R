## Generate loop counts/ call loops for control & NaCl-treated T47D cells processed with `-isDroso true` and `-isDroso false` parameters

## install github
# remotes::install_github("EricSDavis/hictoolsr")

# arsenal -----------------------------------------------------------------

library(hictoolsr)
library(glue)
library(dbscan)
library(tidyverse)
library(InteractionSet)

# loop list set-up --------------------------------------------------------

## vectorize conditions for glue
cond <- c("cont", "sorb", "omega")

## list file paths for `isDroso true` loops
isDroso_loops <- list.files(glue("data/raw/microc/hg38_220801/sip-loops/isDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

## list file paths for `isDroso false` loops
noDroso_loops <- list.files(glue("data/raw/microc/hg38_220801/sip-loops/noDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

## combine them 
loop_files <- c(isDroso_loops, noDroso_loops)

# for extracting counts -----------------------------
hic_files <- list.files(glue("data/raw/microc/hg38_220801/220716_dietJuicerCore/{cond}"),
                        full.names = T)

# merge & extract ---------------------------------------------------------

loopCounts <- 
  mergeBedpe(bedpeFiles = loop_files,
             res = 10000,
             selectCol = 12,
             dist_method = "manhattan",
             minPts = 2) |> 
  as_ginteractions() |> 
  binBedpe(res = 10e3, 
           a1Pos = "center", 
           a2Pos = "center") |> 
  extractCounts(hic = hic_files,
                chroms = paste0("chr",c(1:22, "X", "Y")),
                res = 10e3,
                norm = "NONE",
                matrix = "observed")

# save data -----------------------------
saveRDS(loopCounts, file = "data/processed/microc/YAPP_microc_loopCounts.rds")
