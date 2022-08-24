# Prep for HOMER ----------------------------------------------------------
library(monaLisa)
library(JASPAR2020)
library(TFBSTools)

## prep genome
txdb_human.genome <- makeTxDbFromEnsembl(organism="Homo sapiens",
                                         release=NA,
                                         circ_seqs=NULL,
                                         server="ensembldb.ensembl.org",
                                         username="anonymous", password=NULL, port=0L,
                                         tx_attrib=NULL)
seqlevels(txdb_human.genome)
## Known motifs representing TF binding site preferences
## Extract al vertebrate motifs from JASPAR2020 package as positional weight matrices (PWMs)
pwms <- getMatrixSet(JASPAR2020,
                     opts = list(matrixtype = "PWM",
                                 tax_group = "vertebrates"))

## run HOMER

motif_enrich_HOMER <- calcBinnedMotifEnrR(gr = focal_peaks,
                                          background = filtered_background,
                                          genomedir = txdb_human.genome,
                                          outdir = "data/processed/atac/",
                                          homerfile = findHomer(),
                                          motifFile = pwms,
                                          verbose = TRUE)


