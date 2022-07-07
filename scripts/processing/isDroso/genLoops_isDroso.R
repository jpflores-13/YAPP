## Generate loop counts/ call loops for cont/sorb HEK293 cells processed with `-isDroso` parameter

## install github
# remotes::install_github("EricSDavis/hictoolsr")

## load package
library(hictoolsr)
library(glue)
library(dbscan)

# for merging loops -----------------------------
cond <- c("cont", "sorb", "omega")
loop_files <- list.files(path = glue("data/raw/hic/hg38/sip-loops/isDroso/{cond}"), full.names = T, pattern = "5kbLoops")

# for extracting counts -----------------------------
hic_files <- list.files(path = glue("data/raw/hic/hg38/220627_dietJuicerCore/{cond}"), full.names = T)

# merge / extract -----------------------------
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
              chroms = c(1:22, "X", "Y"),
              res = 10e3,
              norm = "KR",
              matrix = "observed")

# save data -----------------------------
save(loopCounts, file = "data/processed/hic/isDroso/loopCounts.rda")
