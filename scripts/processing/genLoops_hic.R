## Generate loop counts/ call loops for cont/sorb HEK293 cells processed with `-isDroso` parameter

## install github
# remotes::install_github("EricSDavis/hictoolsr")

## load package
library(hictoolsr)
library(glue)
library(dbscan)
library(tidyverse)
library(InteractionSet)

# for merging loops -----------------------------
cond <- c("cont", "sorb", "omega")
isDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/isDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

noDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/noDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

bothDroso_loops <- c(isDroso_loops, noDroso_loops)

# for extracting counts -----------------------------
hic_files <- list.files(path = glue("data/raw/hic/hg38/220722_dietJuicerCore/{cond}"), full.names = T)

# merge / extract -----------------------------
mergedLoops <- 
  mergeBedpe(bedpeFiles = bothDroso_loops,
             res = 10000) |> 
  as_ginteractions()

seqlevelsStyle(mergedLoops) <- "UCSC"

mergedLoops <- mergedLoops |> 
  binBedpe(res = 10e3,
           a1Pos = "center",
           a2Pos = "center")

loopCounts <- extractCounts(bedpe = mergedLoops,
                hic = hic_files,
                chroms = paste0("chr",c(1:22, "X", "Y")),
                res = 10e3,
                norm = "NONE",
                matrix = "observed")
# save data -----------------------------
saveRDS(loopCounts, file = "data/processed/hic/YAPP_hic_loopCounts.rds")

