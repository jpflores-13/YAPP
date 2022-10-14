## Generate loop counts/ call loops for cont/sorb HEK293 cells processed with both `-isDroso` parameters (T/F)

## install github
# remotes::install_github("EricSDavis/hictoolsr")

library(hictoolsr)
library(glue)
library(dbscan)
library(GenomeInfoDb)

# Compile loop files of interest ------------------------------------------

cond <- c("cont", "sorb", "omega")
isDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/isDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

noDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/noDroso/{cond}"),
                            full.names = T,
                            pattern = "5kbLoops")

bothDroso_loops <- c(isDroso_loops, noDroso_loops)
saveRDS(bothDroso_loops, "data/processed/hic/combined_cond_bothDroso_loops.rds")

# Provide .hic files to extract from --------------------------------------

hic_files <- list.files(path = glue("data/raw/hic/hg38/220722_dietJuicerCore/{cond}"), full.names = T)

# Merge bedpe files into one and extract counts  --------------------------

## merge & convert to GInteractions
mergedLoops <- 
  mergeBedpe(bedpeFiles = bothDroso_loops,
             res = 10e3) |> 
  as_ginteractions()

## rename seqstyles level
seqlevelsStyle(mergedLoops) <- "UCSC"

## Bin bedpe files to correct resolution
mergedLoops <- mergedLoops |> 
  binBedpe(res = 10e3,
           a1Pos = "center",
           a2Pos = "center")

loopCounts <- extractCounts(bedpe = mergedLoops,
                hic = hic_files,
                chroms = paste0("chr",c(1:22, "X", "Y")),
                res = 10e3,
                norm = "NONE", ## should always be set to "NONE"
                matrix = "observed")

# save data -----------------------------
saveRDS(loopCounts, file = "data/processed/hic/YAPP_hic_loopCounts.rds")


## save sessionInfo()
sessionInfo()

