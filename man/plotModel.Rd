% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvTools.R
\name{plotModel}
\alias{plotModel}
\title{plotModel}
\usage{
plotModel(
  spec,
  model,
  exclude_zones = NULL,
  labels = c("ppm", "id"),
  groups = NULL,
  tags = FALSE,
  xlab = "",
  ylab = "",
  title = "",
  grp_colors = NULL
)
}
\arguments{
\item{spec}{a 'spec' object (see \code{readSpectrum}, \code{Spec1rDoProc})}

\item{model}{a 'model' object (see \code{specDeconv}, \code{peakOptimize}, \code{GSDeconv}, \code{LSDeconv})}

\item{exclude_zones}{a list of vector defining the excluded zones for lorentzian plots}

\item{labels}{choose as legend labels either 'ppm' or 'id'}

\item{tags}{boolean allowing you to put identification tags at the top of each peak}

\item{title}{title of the graphic}

\item{grp_colors}{specifies the colors for the first groups and/or peaks}
}
\description{
\code{plotModel} plots the model along with the resulting voigt functions from deconvolution
}
