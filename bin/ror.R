#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

tsv <- args[1]
source <- "/fs1/sima/rnaseq_scanb_dev/bin/sourcefiles"
source(paste(source, "applySSP_v1.2.R", sep="/"))

ssp <- "/fs1/sima/rnaseq_scanb_dev/bin/testdata/Training_Run19081Genes_noNorm_SSP.scaled.ROR.tot.asT0.c005.Fcc15_5x5foldCV.num.rules.50_21.selRules.AIMS.GS.RData"
myresults <- applySSP(tsv, ssp, source)
print(myresults)
