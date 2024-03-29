% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doProcessing.R
\name{doProcessing}
\alias{doProcessing}
\title{doProcessing}
\usage{
doProcessing(
  path,
  cmdfile,
  samplefile = NULL,
  bucketfile = NULL,
  phcfile = NULL,
  ncpu = 1
)
}
\arguments{
\item{path}{The full path of either the raw spectra directory on the disk}

\item{cmdfile}{The full path name of the Macro-commands file for processing (text format)}

\item{samplefile}{The full path name of the Sample file (tabular format)}

\item{bucketfile}{The full path name of the file of bucket's zones (tabular format)}

\item{phcfile}{The full path name of the phasing file for samples if required (tabular format)}

\item{ncpu}{The number of cores [default: 1]}
}
\value{
\code{doProcessing} returns a list containing the following components:
\itemize{
  \item \code{samples} : the samples matrix with the correspondence of the raw spectra, 
as well as the levels of the experimental factors if specified in the input.
  \item \code{factors} : the factors matrix with the corresponding factor names. 
At minimum, the list contains the Samplecode label corresponding to the samples without their 
group level.
  \item \code{rawids} : list of the full directories of the raw spectra (i.e. where the FID files 
are accessible)
  \item \code{infos} : list of the acquisition and processing parameters for each (raw) spectra.
  \item \code{specMat} : objects list  regarding the spectra data.
      \itemize{
            \item \code{int} : the matrix of the spectra data (\code{nspec} rows X \code{size} 
columns)
            \item \code{nspec} : the number of spectra
            \item \code{size} : the size (i.e number of points) of each spectra
            \item \code{ppm_min}, \code{ppm_max} : the minimum and the maximum ppm values of 
spectra
            \item \code{ppm} : the vector of the ppm values (\code{size} values)
            \item \code{dppm} : the ppm increment between each point
            \item \code{buckets_zones} : the matrix of the buckets zones including two columns 
(min and max) 
            \item \code{namesASintMax} : boolean - If TRUE, generate all output matrix with 
bucket names based on ppm values of the maximum of the average intensity of all spectra within
the ppm range of each bucket. If FALSE (default), then bucket names will be based on the ppm 
range center of each bucket.
        }
}
}
\description{
\code{doProcessing} is the main function of this package. Indeed, this function performs 
the complete processing of a set of 1D NMR spectra from the FID (raw data) and based on a 
processing sequence (macro-command file). An additional file specifies all the spectra to 
be considered by associating their sample code as well as the levels of experimental 
factors to which they belong. In this way it is possible to select only a subset of spectra 
instead of the whole set.
}
\examples{
 \donttest{
    data_dir <- system.file("extra", package = "Rnmr1D")
    cmdfile <- file.path(data_dir, "NP_macro_cmd.txt")
    samplefile <- file.path(data_dir, "Samples.txt")
    out <- Rnmr1D::doProcessing(data_dir, cmdfile=cmdfile, 
                                samplefile=samplefile, ncpu=2)
}
}
\seealso{
the NMRProcFlow online documentation \url{https://nmrprocflow.org/} and especially 
the Macro-command Reference Guide (\url{https://nmrprocflow.org/themes/pdf/Macrocommand.pdf})
}
