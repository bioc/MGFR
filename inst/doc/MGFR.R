### R code from vignette source 'MGFR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: MGFR.Rnw:56-58
###################################################
library("MGFR")
ls("package:MGFR")


###################################################
### code chunk number 2: MGFR.Rnw:112-113
###################################################
options(width=60)


###################################################
### code chunk number 3: MGFR.Rnw:115-118
###################################################
data("ref.mat")
dim(ref.mat)
colnames(ref.mat)


###################################################
### code chunk number 4: reference
###################################################
library(MGFR)
data("ref.mat")
markers.list <- getMarkerGenes.rnaseq(ref.mat, class.vec = colnames(ref.mat), samples2compare="all", annotate=TRUE, gene.ids.type="ensembl", score.cutoff=1)
names(markers.list)
# show the first 20 markers of liver
markers.list[["liver_markers"]][1:20]


###################################################
### code chunk number 5: sessionInfo
###################################################
sessionInfo()


