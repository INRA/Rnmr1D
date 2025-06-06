% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvTools.R
\name{getSlices}
\alias{getSlices}
\title{getSlices}
\usage{
getSlices(spec, ppmranges, flatwidth = 0.004, snrfactor = 4, maxrange = 0.3)
}
\arguments{
\item{spec}{a 'spec' object}

\item{ppmranges}{ppm ranges as a matrix in order to apply the deconvolution, each row specifying a zone}

\item{flatwidth}{specifies the minimum width of a zone in which the spectrum intensities are close to zero to consider this one as a cutting zone (default 0.003 ppm)}

\item{snrfactor}{specifies factor applied on the Std Dev. of the noise used as threshold for first derivate intensities (default=4)}

\item{maxrange}{specifies the maximum width of a cutting zone (default 0.3 ppm)}
}
\value{
a list of ppm range
}
\description{
Slice the spectrum in order to define ranges for Local Deconvolution (LSDeconv) and return only those include the provided ppmranges
}
