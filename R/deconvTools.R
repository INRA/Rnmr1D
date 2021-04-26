# ID deconvTools.R
# Copyright (C) 2017-2021 INRAE
# Authors: D. Jacob



#' filterSavGol
#'
#' \code{filterSavGol} applies a Savitzky-Golay filter on a spectral signal.
#' @param s the spectral signal as a numerical vector
#' @param m the degree of the polynomial filter (integer)
#' @param nl width of the sliding window on the left (integer)
#' @param nr width of the sliding window on the rigth (integer)
#' @return a vector of the same dimension as the entry one
filterSavGol <- function(s, m, nl, nr)
{
   C_fSavGol(s, m, nl, nr)
}


#' filterByWT
#'
#' \code{filterByWT} applies a filtering based on wavelet using the universal threshold
#' @param s the spectral signal as a numerical vector
#' @param wavelet the type of the wavelet: 0 for HAAR, 1 for DAUB2, 2 for DAUB4, 3 for DAUB8
#' @param type the type of the threshold :  0 for Soft threshold, 1 for Hard threshold
#' @return a vector of the same dimension as the entry one
filterByWT <- function(s, wavelet, type = 0)
{
   C_FilterbyThreshold(s, wavelet, threshold = type)
}


#' filterByThreshold
#'
#' \code{filterByThreshold} applies a filtering based on wavelet by specifying a threshold value 
#' @param s the spectral signal as a numerical vector
#' @param wavelet the type of the wavelet: 0 for HAAR, 1 for DAUB2, 2 for DAUB4, 3 for DAUB8
#' @param threshold the threshold value - default value is 0.5
#' @return a vector of the same dimension as the entry one
filterByThreshold <- function(s, wavelet, threshold = 0.5)
{
   C_FilterbyWT(s, type = wavelet, threshold)
}


#' Lorentz
#'
#' \code{Lorentz} belongs to the low-level functions group for deconvolution.
#' @param ppm a vector of ppm values
#' @param amp amplitude of the lorentzian
#' @param x0 central value of the lorentzian
#' @param sigma half-width of the lorentzian
#' @return a vector of the lorentzian values (same size as ppm)
Lorentz <- function(ppm, amp, x0, sigma)
{
   C_Lorentz(ppm, amp, x0, sigma)
}


#' optimOneLorentz
#'
#' \code{optimOneLorentz} belongs to the low-level functions group for deconvolution.
#' @param X a vector of ppm values
#' @param Y a vector of intensities
#' @param par a vector of the 3 lorentzian parameters namely: Amplitude, central ppm value, ppm width at mid-height
#' @return a vector of the lorentzian parameters (same size as par)
optimOneLorentz <- function(X, Y, par)
{
   C_OneLorentz(X, Y, par)
}


#' MyFuncTest
#'
#' \code{MyFuncTest} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution 
#' @param filt a list of specific parameters for filtering
#' @param peaks a list of specific parameters for deconvolution
#' @param verbose level of debug information
#' @return a list
MyFuncTest <- function(spec, ppmrange, filt, peaks, verbose=1)
{
   C_MyFuncTest(spec, ppmrange, filt, peaks, verbose)
}


#' MyFuncTest2
#'
#' \code{MyFuncTest2} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param n1 start index of the spectrum segment to be corrected
#' @param n2 stop index of the spectrum segment to be corrected
#' @return a vector of the same dimension as spec$int
MyFuncTest2 <- function(spec, n1, n2)
{
   C_MyFuncTest2(spec, n1, n2)
}


