#' Lorentz
#'
#' \code{Lorentz} belongs to the low-level functions group for deconvolution.
#' @param ppm a vector of ppm values
#' @param amp amplitude of the lorentzian
#' @param x0 central value of the lorentzian
#' @param sigma half-width of the lorentzian
#' @return a vector of the lorentzian values (same size as ppm)
Lorentz = function(ppm, amp, x0, sigma)
{
   C_Lorentz(ppm, amp, x0, sigma)
}

#' MyFuncTest
#'
#' \code{MyFuncTest} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list of two values
#' @param filt a list of specific parameters for filtering
#' @param peaks a list of specific parameters for deconvolution
#' @param verbose level of debug information
#' @return a list
MyFuncTest = function(spec, ppmrange, filt, peaks, verbose=1)
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
MyFuncTest2 = function(spec, n1, n2)
{
   C_MyFuncTest2(spec, n1, n2)
}
