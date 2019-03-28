\name{Rnmr1D-package}
\alias{Rnmr1D-package}
\alias{Rnmr1D}
\docType{package}
\title{
  Rnmr1D is the main module in the NMRProcFlow web application (nmrprocflow.org) concerning the NMR spectra processing.
}
\description{
Rnmr1D R package is aimed to performs the complete processing of a set of 1D NMR spectra from the FID (raw data) and based on a processing sequence (macro-command file). An additional file specifies all the spectra to be considered by associating their sample code as well as the levels of experimental factors to which they belong.
}
\details{
  NMRProcFlow allows experts to build their own spectra processing workflow, in order to become re-applicable to similar NMR spectra sets, i.e. stated as use-cases. By extension, the implementation of NMR spectra processing workflows executed in batch mode can be considered as relevant provided that we want to process in this way very well-mastered and very reproducible use cases, i.e. by applying the same Standard Operating Procedures (SOP). A subset of NMR spectra is firstly processed in interactive mode in order to build a well-suited workflow. This mode can be considered as the 'expert mode'. Then, other subsets that are regarded as either similar or being included in the same case study, can be processed in batch mode, operating directly within a R session.
  
  
  This package consists of two sets of functions depending on the use: a so-called "high level" set and a so-called "low level" set.
  
  * The "high level" group consists of the following functions :

       \code{\link{detectCores}}, \code{\link{doProcessing}}, \code{\link{doProcCmd}}, \code{\link{getBucketsDataset}}, \code{\link{getBucketsTable}}, \code{\link{getClusters}}, \code{\link{getMergedDataset}}, \code{\link{getSnrDataset}}, \code{\link{getSpectraData}}, \code{\link{plotClusters}}, \code{\link{plotCriterion}}, \code{\link{plotLoadings}},  \code{\link{plotScores}}, \code{\link{plotSpecMat}}

    to which must be added the graphical output methods based on ggplot2 and plotly :

       \code{\link{ggplotClusters}}, \code{\link{ggplotCriterion}}, \code{\link{ggplotLoadings}},  \code{\link{ggplotScores}}, \code{\link{ggplotPlotly}}
    
    This set of function, dedicated to 1D NMR spectra  (1H & 13C) for metabolomics  allow user to process  a complet set of 1D NMR spectra from the FID (raw data) and visualize the resulting spectra in the frequency domain. In output, a plenty of data are available : Bucket table, Data matrix, Signal-to-noise ratios matrix,  spectral data, metadata on both samples and factors. It is possible to reduce the variables size by gathering within clusters two or more buckets  likely to be the (partial) spectral signature of a single chemical compound. The example below provides an illustration of a typical use.

  * The "low level" group consists of the following functions :
 
       \code{\link{generateMetadata}}, \code{\link{Spec1rProcpar}}, \code{\link{Spec1rDoProc}}, \code{\link{RWrapperCMD1D}}

    This set of functions allows to process spectra one by one (Spec1rDoProc) and to apply the treatments one by one on the spectra set ('RWrapperCMD1D'). Being intended to be integrated into a complex workflow such as implemented in the NMRProcFlow application, these functions are accessible only to maintain compatibility between this package and the NMRProcFlow application, and are therefore very little documented.
}

\author{
Daniel Jacob

Maintainer: Daniel Jacob <daniel.jacob@inra.fr>
}
\references{
   Jacob, D., Deborde, C., Lefebvre, M., Maucourt, M. and Moing, A. (2017) NMRProcFlow: A graphical and interactive tool dedicated to 1D spectra processing for NMR-based metabolomics, Metabolomics 13:36. \url{https://www.ncbi.nlm.nih.gov/pubmed/28261014}
}
\keyword{ package }
\seealso{
  the NMRProcFlow online documentation \url{https://nmrprocflow.org/}
}
\examples{
  \dontrun{
     ## Process the raw data 
     data_dir <- system.file("extra", package = "Rnmr1D")
     RAWDIR <- file.path(data_dir, "CD_BBI_16P02")
     CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
     SAMPLEFILE <- file.path(data_dir, "Samples.txt")
     out <- Rnmr1D::doProcessing(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=detectCores())

     ## Plot the spectra for the given ppm range
     plotSpecMat(out$specMat, ppm_lim=c(0.5,9), K=0)

     ## Generate the output data matrix based on the computed buckets 
     ## stored in  out$specMat$bucket_zones.
     outMat <- Rnmr1D::getBucketsDataset(out, norm_meth='CSN')

     ## Compute the bucket clustering based on a cut threshold  applied on the buckets tree 
     ## resulting of a hierarchical classification analysis. If \code{vcutusr} is equal to 
     ## zero then cut threshold automatically estimated; otherwise the specified value 
     ## is applied.
     outclust <- Rnmr1D::getClusters(outMat, method='hca', vcutusr=0.15 )

     ## Plot the criterion curves and the cluster boxplot
     dev.new()
     plotCriterion(outclust)
     dev.new()
     plotClusters(outMat,outclust)

     ## ---- PCA  ----
     pca <- prcomp(outMat,retx=TRUE,scale=TRUE, rank=2)

     # plot PCA Scores
     dev.new()
     plotScores(pca$x, 1, 2, out$samples, 'Treatment', level=0.95)  # Choose 'Treatment' as factor

     # plot PCA Loadings
     dev.new()
     plotLoadings(pca$rotation, 1, 2, associations=outclust$clustertab, cexlabel=0.6, 
                  main=sprintf("Loadings - Crit=%s",outclust$vcrit) )

  }
}
