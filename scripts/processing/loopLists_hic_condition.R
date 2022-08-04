library(hictoolsr)
library(InteractionSet)
library(tidyverse)
library(glue)
library(dbscan)
library(data.table)

# create loop lists --------------------------------------------------------
## Merge gained control loops from both `-isDroso true` & `-isDroso false`
cont_isDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/isDroso/cont"),
                                 full.names = T,
                                 pattern = "5kbLoops")

cont_noDroso_loops <- list.files(path = glue("data/raw/hic/hg38/sip-loops/noDroso/cont"),
                                 full.names = T,
                                 pattern = "5kbLoops")

## Saved as RDS to speed-up analysis
cont_loops <- c(cont_noDroso_loops, cont_isDroso_loops) |>
  mergeBedpe(res = 10e3,
             selectCol = 12,
             dist_method = "manhattan",
             minPts = 2)
saveRDS(cont_loops, "data/processed/hic/cont_bothDroso_loops.rds")

## Merge sorbitol loops from both `-isDroso true` & `-isDroso false`
sorb_isDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/isDroso/sorb",
                                 full.names = T,
                                 pattern = "5kbLoops")

sorb_noDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/noDroso/sorb",
                                 full.names = T,
                                 pattern = "5kbLoops")

## Saved as RDS to speed-up analysis
sorb_loops <- c(sorb_isDroso_loops, sorb_noDroso_loops) |>
  mergeBedpe(res = 10e3,
             selectCol = 12,
             dist_method = "manhattan",
             minPts = 2)
saveRDS(sorb_loops, "data/processed/hic/sorb_bothDroso_loops.rds")

## Merge omega loops form both `-isDroso true` & `-isDroso false`
omega_isDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/isDroso/omega/",
                                  full.names = T,
                                  pattern = "5kbLoops")

omega_noDroso_loops <- list.files(path = "data/raw/hic/hg38/sip-loops/noDroso/omega/",
                                  full.names = T,
                                  pattern = "5kbLoops")

## Saved as RDS to speed-up analysis
omega_loops <- c(omega_isDroso_loops, omega_noDroso_loops) |>
  mergeBedpe(res = 10e3,
             selectCol = 12,
             dist_method = "manhattan",
             minPts = 2)
saveRDS(omega_loops, "data/processed/hic/omega_bothDroso_loops.rds")
