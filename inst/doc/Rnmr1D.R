## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----init-object0, echo=TRUE, include=TRUE-------------------------------
library(Rnmr1D)
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "MMBBI_14P05")
CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
SAMPLEFILE <- file.path(data_dir, "Samples.txt")

## ----init-object1, echo=TRUE, include=TRUE-------------------------------
samples <- read.table(SAMPLEFILE, sep="\t", header=T,stringsAsFactors=FALSE)
samples

## ----init-object2, echo=TRUE, include=TRUE------------------------------------------------------------------------------------
CMDTEXT <- readLines(CMDFILE)
options(width=128)
CMDTEXT[grep("^#$", CMDTEXT, invert=TRUE)]

## ----eval0, echo=TRUE, eval=FALSE---------------------------------------------------------------------------------------------
#  out <- Rnmr1D::doProcessing(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=detectCores())

## ----proc1, echo=FALSE, eval=TRUE---------------------------------------------------------------------------------------------
data(test)
out <- test[[1]]
outclust <- test[[2]]

## ----proc2, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
ls(out)

out$infos

## ----plot1, echo=TRUE, fig.align='center', fig.width=12-----------------------------------------------------------------------
plotSpecMat(out$specMat, ppm_lim=c(0.5,5))

## ----plot2, echo=TRUE, fig.align='center', fig.width=12-----------------------------------------------------------------------
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0)

## ----plot3, echo=TRUE, fig.align='center', fig.width=12-----------------------------------------------------------------------
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33, asym=0)
cols <- c(rep("red",6), rep("blue",6))
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.67, dppm_max=0, cols=cols)

## ----proc3, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
outMat <- Rnmr1D::getBucketsDataset(out, norm_meth='CSN')
outMat[, 1:10]

## ----proc4, echo=TRUE, eval=FALSE---------------------------------------------------------------------------------------------
#  outclust <- Rnmr1D::getClusters(outMat, method='corr', cval=0, dC=0.01, ncpu=detectCores())

## ----proc4b, echo=TRUE, eval=TRUE---------------------------------------------------------------------------------------------
outclust$clustertab[1:20, ]
outclust$clusters$C4      # same as outclust$clusters[['C4']]

## ----plot4, echo=TRUE, fig.align='center', fig.width=12, fig.height=8---------------------------------------------------------
plotCriterion(outclust)

## ----plot5, echo=TRUE, fig.align='center', fig.width=12, fig.height=8---------------------------------------------------------
plotClusters(outMat,outclust)

## ----proc5, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
pca <- prcomp(outMat,retx=TRUE,scale=T, rank=2)

## ----plot6, echo=TRUE, fig.align='center', fig.width=12, fig.height=8---------------------------------------------------------
plotScores(pca$x, out$samples, 'Genotype', level=0.8)   # Choose 'Genotype' as factor, confidence level = 80%

plotScores(pca$x, out$samples, 'Treatment', level=0.8)  # Choose 'Treatment' as factor, confidence level = 80%

## ----plot7, echo=TRUE, fig.align='center', fig.width=12, fig.height=10--------------------------------------------------------
plotLoadings(pca$rotation, 1, 2, associations=outclust$clustertab, cexlabel=0.6, main=sprintf("Loadings - Crit=%s",outclust$vcrit) )

## ----plot8, echo=TRUE, fig.align='center', fig.width=12, fig.height=10--------------------------------------------------------
outMat.merged <- Rnmr1D::getMergedDataset(outMat, outclust)
pca.merged <- prcomp(outMat.merged,retx=TRUE,scale=T, rank=2)
plotLoadings(pca.merged$rotation, 1,2,  associations=outclust$clustertab, cexlabel=0.6, main=sprintf("Loadings - Crit=%s",outclust$vcrit) )

