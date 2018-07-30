## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----init-object0, echo=TRUE, include=TRUE-------------------------------
library(Rnmr1D)
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "CD_BBI_16P02")
CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
SAMPLEFILE <- file.path(data_dir, "Samples.txt")

## ----init-object1, echo=TRUE, include=TRUE-------------------------------
samples <- read.table(SAMPLEFILE, sep="\t", header=T,stringsAsFactors=FALSE)
samples

## ----init-object2, echo=TRUE, include=TRUE------------------------------------------------------------------------------------
CMDTEXT <- readLines(CMDFILE)
options(width=128)
CMDTEXT[grep("^#$", CMDTEXT, invert=TRUE)]

## ----eval0, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
out <- Rnmr1D::doProcessing(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=2)

## ----proc2, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
ls(out)

out$infos

## ----plot1, echo=TRUE, fig.align='center', fig.width=12-----------------------------------------------------------------------
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33)

## ----plot2, echo=TRUE, fig.align='center', fig.width=12-----------------------------------------------------------------------
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0, pY=0.1)

## ----plot3, echo=TRUE, fig.align='center', fig.width=12-----------------------------------------------------------------------
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33, asym=0)
cols<- c( rep("blue",length(out$samples$Treatment)));
cols[out$samples$Treatment=="stress"] <- "red"
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.67, dppm_max=0, cols=cols)

## ----proc3, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
specMat.new <- Rnmr1D::doProcCmd(out, 
     c( "bucket aibin 10.2 10.5 0.3 3 0", "9.5 4.9", "4.8 0.5", "EOL" ), ncpu=2, debug=TRUE)
out$specMat <- specMat.new

## ----proc4, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
outMat <- Rnmr1D::getBucketsDataset(out, norm_meth='CSN')
outMat[, 1:10]

## ----proc5, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
options(warn=-1)
outclust <- Rnmr1D::getClusters(outMat, method='corr', cval=0, dC=0.005, ncpu=2)

## ----proc6, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
outclust$clustertab[1:20, ]
outclust$clusters$C5      # same as outclust$clusters[['C5']]

## ----plot5, echo=TRUE, fig.align='center', fig.width=12, fig.height=8---------------------------------------------------------
plotCriterion(outclust)

## ----plot6, echo=TRUE, fig.align='center', fig.width=12, fig.height=8---------------------------------------------------------
plotClusters(outMat,outclust)

## ----proc7, echo=TRUE, eval=TRUE----------------------------------------------------------------------------------------------
pca <- prcomp(outMat,retx=TRUE,scale=T, rank=2)

## ----plot7, echo=TRUE, fig.align='center', fig.width=12, fig.height=8---------------------------------------------------------
plotScores(pca$x, 1, 2, out$samples, factor='Treatment', level=0.95)  # Choose 'Treatment' as factor, confidence level = 95%

## ----plot8, echo=TRUE, fig.align='center', fig.width=12, fig.height=10--------------------------------------------------------
plotLoadings(pca$rotation, 1, 2, associations=outclust$clustertab, 
             cexlabel=0.6, level=0.6, main=sprintf("Loadings - Crit=%s",outclust$vcrit) )

## ----plot9, echo=TRUE, fig.align='center', fig.width=12, fig.height=10--------------------------------------------------------
outMat.merged <- Rnmr1D::getMergedDataset(outMat, outclust)
pca.merged <- prcomp(outMat.merged,retx=TRUE,scale=T, rank=2)
plotLoadings(pca.merged$rotation, 1, 2, associations=outclust$clustertab, 
             cexlabel=0.6, main=sprintf("Loadings - Crit=%s",outclust$vcrit) )

## ----proc101, echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "CD_BBI_16P02")

## ----proc102, echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------
procParams <- Spec1rProcpar
procParams$LOGFILE <- ""
procParams$VENDOR <- 'bruker'
procParams$INPUT_SIGNAL <- 'fid'
procParams$LB <- 0.3
procParams$ZEROFILLING <- TRUE
procParams$ZFFAC <- 4
procParams$OPTPHC1 <- TRUE
procParams$TSP <- TRUE

## ----proc103, echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------
metadata <- generateMetadata(RAWDIR, procParams)
metadata

## ----proc104, echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------
ID <- 1
ACQDIR <- metadata$rawids[ID,1]
spec <- Spec1rDoProc(Input=ACQDIR,param=procParams)

## ----proc105, echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------
ls(spec)

## ----plot100, echo=TRUE, fig.align='center', fig.width=12, fig.height=8-------------------------------------------------------
plot( spec$ppm, spec$int, type="l", col="blue", 
                xlab="ppm", ylab="Intensities", 
                xlim=c( spec$pmax, spec$pmin ), ylim=c(0, max(spec$int/100)) )
legend("topleft", legend=metadata$samples[ID,1])

