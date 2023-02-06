# ID deconvTools.R
# Copyright (C) 2017-2022 INRAE
# Authors: D. Jacob
#
#=====================================================================
# Deconvolution parameters
#=====================================================================


# Filter types
fnone <- list( type=0 )
fdaub8 <- list( type=1, threshold=0.5 )
fsymlet8 <- list( type=2, threshold=0.5 )
fsavgol3 <- list( type=3, m=3, nl=5, nr=5 )
fsavgol5 <- list( type=3, m=5, nl=8, nr=8 )
fsavgol10 <- list( type=3, m=10, nl=20, nr=30 )
fsavgol50 <- list( type=3, m=50, nl=60, nr=80 )
fsmooth <- list( type=4, m=12 )

# Names of Filter type
filtnames <- list('haar'=0, 'daub2'=1, 'daub4'=2, 'daub8'=3, 'symlet2'=4, 'symlet4'=5, 'symlet8'=6)

#' deconvParams
#'
#' Initialize the deconvolution parameter list
#' @return
#' \itemize{
#'   \item \code{flist} : Filter type list : 'smooth1', 'smooth2' and'smooth3' for Savitzky-Golay filter, 'daub8' and 'symlet8' for filter based on wavelet
#'   \item \code{criterion} : Criterion type for the optimizations : 0 => R2, 1 => 1/Std(residues) - default value = 0
#'   \item \code{reltol} : Criterion tolerance for the optimization - default value = 0.0001
#'   \item \code{facN} :  Noise factor applied while the peak finding step - default value = NULL
#'   \item \code{ratioPN} : Peak/Noise Ratio applied while the peak selection step - default value = 1
#'   \item \code{obl} : Optimization of a baseline (BL) for each massif. 0 means no BL, an integer greater than 0 indicates the polynomial order of the BL default value = 0
#'   \item \code{distPeaks} : PeakFinder : min distance between 2 peaks (as multiple of sigma_min which is typically equal to 0.0005 ppm) - default value = 2
#'   \item \code{optim} : Indicates if optimisation is applied - default value = 1
#'   \item \code{oppm} : Indicates if ppm optimisation is applied - default value = 1
#'   \item \code{osigma} : Indicates if sigma optimisation is applied - default value = 1
#'   \item \code{d2meth} : PeakFinder : Indicates if minima method to the second derivation is applied
#'   \item \code{spcv} : PeakFinder : Maximal CV on Spectrum - default value = 0.005
#'   \item \code{d2cv} : PeakFinder : Maximum CV on the derived spectrum - default value = 0.05
#'   \item \code{d1filt} : Apply Filter (1) on the 1st derivate or not (0) - default value = 0
#'   \item \code{d2filt} : Apply Filter (1) on the 2nd derivate or not (0) - default value = 0
#'   \item \code{sigma_min} : Optimization of Sigmas : Fixe the minimum limit of sigmas - default value = 0.0005
#'   \item \code{sigma_max} : Optimization of Sigmas : Fixe the maximum limit of sigmas - default value = 0.005
#'   \item \code{verbose} : Indicates if we want information messages - default value = 1
#'   \item \code{exclude_zones} : Exclude ppm zones for the criterion evaluation - default value = NULL
#' }
deconvParams <- list (
  # Filter types
  flist = list( 'smooth0'=fsavgol3, 'smooth1'=fsavgol5, 'smooth2'=fsavgol10, 'smooth3'=fsavgol50,
                'daub8'=fdaub8, 'symlet8'=fsymlet8, 'none'=fnone ),

  # Criterion type for the optimizations
  # 0 => R2
  # 1 => 1/Std(residues)
  criterion = 0,

  # Criterion tolerance for the optimization
  reltol = 0.0001,
  maxstep = 50,

  # Peak/Noise Ratio
  facN = NULL,
  ratioPN = 5,
  lowPeaks = 1,
  
  # indicates if pseudo-voigt is used instead of lorentzian
  pvoigt=1,
  eta=0.7,
  asym=0,
  asymmax=50,

  # Optimization of peaks : 0 => No, 1 => Yes
  optim=1,
  oppm = 1,
  osigma = 1,
  oasym = 0,
  oeta=1,
  estimate_int=0,

  # Indicate if peaks are selected based on the threshold ratio Signal-Noise
  selectpk = 0,

  # Optimization by only one block or by several blocks applying a cut-off process. 
  oneblk=1,
  scmin=2,

  # Optimization of a baseline (BL) for each massif
  # 0 means no BL, an integer greater than 0 indicates the polynomial order of the BL
  obl = 0,

  # Peaks searching : min distance between 2 peaks (as multiple of sigma_min which is typically equal to 0.0005 ppm)
  distPeaks = 2,

  # Peaks searching : Minima method applied to the second derivation
  d2meth = 1,

  # minimum values for the peak curvature : the larger the value, the narrower the peak 
  # 1) on spectrum (sp) and 2) its second derivate (d2)
  spcv = 0.01, # 0.005 - 0.025
  d2cv = 0.1,  # 0.05 - 0.25

  # Apply Filter (1) or not (0) on 1st (d1filt) and 2nd derivates (d2filt)
  d1filt = 0,
  d2filt = 0,

  # Optimization of Sigmas : Fixe the limits (min and max) of sigmas
  sigma_min = 0.0005,
  sigma_max = 0.005,

  # Indicates if we want information messages
  verbose = 1,

  # Exclude ppm zones for the criterion evaluation
  exclude_zones = NULL,
  
  # a dataframe of peaks (columns : pos, ppm, amp, sigma, eta, integral)
  peaks = NULL,

  # Indicates if a second deconvolution phase without the highest peaks has to be done
  spass = 0,

  # slicing: the minimum width of a zone in which the spectrum intensities are close to zero to consider this one as a cutting zone
  flatwidth = 0.004,

  # slicing: the factor applied on the Std Dev. of the noise used as threshold for first derivate intensities
  snrfactor = 4,

  # slicing: the maximum width of a cutting zone (default 0.3 ppm)
  maxrange = 0.3
)

#=====================================================================
# Filtering
#=====================================================================

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
#' @param wavelet the name of the wavelet: haar, daub2, daub4, daub8, symlet2, symlet4, symlet8
#' @param type the type of the threshold :  0 for Soft threshold, 1 for Hard threshold
#' @return a vector of the same dimension as the entry one
filterByWT <- function(s, wavelet, type = 0)
{
   C_FilterbyThreshold(s, filtnames[[wavelet]], threshold = type)
}


#' filterByThreshold
#'
#' \code{filterByThreshold} applies a filtering based on wavelet by specifying a threshold value
#' @param s the spectral signal as a numerical vector
#' @param wavelet the name of the wavelet: haar, daub2, daub4, daub8, symlet2, symlet4, symlet8
#' @param threshold the threshold value - default value is 0.5
#' @return a vector of the same dimension as the entry one
filterByThreshold <- function(s, wavelet, threshold = 0.5)
{
   C_FilterbyWT(s, type = filtnames[[wavelet]], threshold)
}

#=====================================================================
# Deconvolution - general routines
#=====================================================================

#' getDeconvParams
#'
#' \code{getDeconvParams} merges some specific parameters values with the full deconvolution list and return the resulting list. With no parameter as input, it returns the default parameter list.
#' @param params a list defining some specific parameters for deconvolution
#' @return the resulting list of deconvolution parameters.
getDeconvParams <- function(params=NULL)
{
  fullparams <- deconvParams
  if ("list" %in% class(params))
     for (p in ls(params)) fullparams[[p]] <- params[[p]]
  fullparams
}

# get the index sequence corresponding to the ppm range
getIndexSeq <- function(spec, ppm)
{
   c(which(spec$ppm>=ppm[1])[1]:length(which(spec$ppm<=ppm[2])))
}

#' pos2ppm
#'
#' \code{pos2ppm} convert an index position to the corresponding ppm value
#' @param spec a 'spec' object
#' @param index an index position
#' @return the corresponding ppm value
pos2ppm <- function(spec, index)
{
   if (index<1 || index>length(spec$ppm))
      stop("the index is out of range")
   return( spec$pmin + (index-1)*spec$dppm )
}

#' ppm2pos
#'
#' \code{ppm2pos} convert a ppm value to the corresponding index position 
#' @param spec a 'spec' object
#' @param ppm a ppm value
#' @return the corresponding index position 
ppm2pos <- function(spec, ppm)
{
   if (ppm<spec$pmin || ppm>spec$pmax)
      stop("the ppm is out of range")
   return( round((ppm - spec$pmin)/spec$dppm,0) + 1)
}

#' Lorentz
#'
#' \code{Lorentz} belongs to the low-level functions group for deconvolution.
#' @param ppm a vector of ppm values
#' @param amp amplitude of the lorentzian
#' @param x0 central value of the lorentzian
#' @param sigma half-width of the lorentzian
#' @param asym asymetric parameter
#' @return a vector of the lorentzian values (same size as ppm)
Lorentz <- function(ppm, amp, x0, sigma, asym)
{
   C_Lorentz(ppm, amp, x0, sigma, asym)
}

#' PVoigt
#'
#' \code{PVoigt} belongs to the low-level functions group for deconvolution.
#' @param ppm a vector of ppm values
#' @param amp amplitude of the lorentzian
#' @param x0 central value of the lorentzian
#' @param sigma half-width of both lorentzian and gaussian
#' @param asym asymetric parameter
#' @param eta mixing coefficient for the pseudo-voigt function (between 0 and 1)
#' @return a vector of the lorentzian values (same size as ppm)
PVoigt <- function(ppm, amp, x0, sigma, asym, eta)
{
   C_PVoigt(ppm, amp, x0, sigma, asym, eta)
}

#' optimOneVoigt
#'
#' \code{optimOneVoigt} belongs to the low-level functions group for deconvolution.
#' @param X a vector of ppm values
#' @param Y a vector of intensities
#' @param par a vector of the 3 pseudo-voigt parameters namely: Amplitude, central ppm value, 2 ppm widths at mid-height for mixed lorentizian and gaussian
#' @return a vector of the pseudo-voigt parameters (same size as par)
optimOneVoigt <- function(X, Y, par)
{
   C_OneVoigt(X, Y, par)
}


#' peakFinder
#'
#' \code{peakFinder} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param params a list of specific parameters for deconvolution
#' @param filter a filter type for filtering the noise and  smoothing the signal
#' @param verbose level of debug information
#' @return a list
peakFinder <- function(spec, ppmrange, params=NULL, filter='none', verbose=1)
{
   g <- getDeconvParams(params)
   model <- C_peakFinder(spec, ppmrange, g$flist[[filter]], g, verbose)
   class(model) = "peakModel"
   model
}

#' peakOptimize
#'
#' \code{peakOptimize} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param params a list of specific parameters for deconvolution, including the matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, eta
#' @param verbose level of debug information
#' @return a list
peakOptimize <- function(spec, ppmrange, params, verbose=1)
{
   if (is.null(params$peaks) || ! "data.frame" %in% class(params$peaks) || nrow(params$peaks)==0 )
      stop("the peaks param must be a data.frame with at least one row")
   model <- C_peakOptimize(spec, ppmrange, params, verbose)
   class(model) = "optimModel"
   model
}

#' peakFiltering
#'
#' \code{peakFiltering} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param peaks the matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, eta
#' @param ratioPN the ratio Peaks/Noise for filtering
#' @return a dataframe
peakFiltering <- function(spec, peaks, ratioPN)
{
   if (is.null(peaks) || ! "data.frame" %in% class(peaks) || nrow(peaks)==0 )
      stop("the peaks param must be a data.frame with at least one row")
   C_peakFiltering(spec, peaks, ratioPN)
}

#' specModel
#'
#' \code{specModel} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param peaks a matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, eta, integral
#' @return a vector with the same size of the spectrum
specModel <- function(spec, ppmrange, peaks)
{
   C_specModel(spec, ppmrange, peaks)
}

#' specDeconv
#'
#' \code{specDeconv} = \code{peakFinder} + \code{peakOptimize}. \code{specDeconv} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param filter a filter type for filtering the noise and  smoothing the signal
#' @param params a list of specific parameters for deconvolution
#' @param verbose level of debug information
#' @return a list
specDeconv <- function(spec, ppmrange, params=NULL, filter='none', verbose=1)
{
   model0 <- peakFinder(spec, ppmrange, params, filter, verbose = verbose)
   params$peaks <- model0$peaks
   peakOptimize(spec, ppmrange, params, verbose = verbose)
}

# Estimated number of peaks
estimation_nbpeaks <- function(spec, ppmrange, params=NULL)
{
   g <- getDeconvParams(params)
   FacN <- ifelse(is.null(g$facN), 5, g$facN)
   spec$B <- spec$Noise/FacN
   g$ratioSN <- ifelse(g$lowPeaks==0, FacN*g$ratioPN, FacN/10)
   model0 <- peakFinder(spec, ppmrange, g, 'daub8', verbose = 0)
   model0$nbpeak
}

#=====================================================================
# Processing
#=====================================================================

# Read then process the spectrum

# directory containing the FID

#' readSpectrum
#'
#' \code{read_spectrum} belongs to the low-level functions group - it processes only one raw spectrum at time.
#' @param ACQDIR Full directory path of the raw spectrum, i.e the directory containing the FID
#' @param procParams a Spec1rProcpar list
#' @param ppmnoise the ppm range containing only noise signal in order to estimate the level of the noise (S/N ratio)
#' @param PHC the phasing values applied to the spectrum in the frequency domain thus avoiding the automatic phasing step. Only useful if the input signal is an FID (procParams$INPUT_SIGNAL)
#' @param scaleIntensity factor applied to the intensities to establish a change of scale.
#' @param verbose level of debug information
#' @return spec object
readSpectrum <- function(ACQDIR, procParams, ppmnoise=c(10.2,10.5), PHC=NULL, scaleIntensity=1, verbose=1)
{
   if (!is.null(PHC)) {
      procParams$OPTPHC0 <<- FALSE
      procParams$OPTPHC1 <<- FALSE
      procParams$phc0 <<- PHC[1]*pi/180
      procParams$phc1 <<- PHC[2]*pi/180
   }

   procParams$DEBUG <- ifelse(verbose, TRUE, FALSE)
   spec <- Rnmr1D::Spec1rDoProc(Input=ACQDIR,param=procParams)
   spec <- Rnmr1D:::.ppm_calibration(spec)
   Noise <- C_noise_estimation(spec$int, which(spec$ppm>=ppmnoise[1])[1], length(which(spec$ppm<=ppmnoise[2])))
   spec$int <- spec$int/scaleIntensity
   spec$B <- spec$Noise <- Noise/scaleIntensity
   spec$size <- length(spec$int)
   if (verbose) {
      cat('Size =',length(spec$int),', Max Intensity =',round(max(spec$int),3),', Noise Level =',round(spec$Noise,6),"\n")
   }
   spec
}

# PPM calibration of the spectrum

#' calibSpectrum
#'
#' \code{calibSpectrum} belongs to the low-level functions group - it processes only one raw spectrum at time.
#' @param spec a spectrum object returned by the readSpectrum function
#' @param zoneref the ppm range containing the TSP/DSS signal
#' @param ppmref the ppm reference value 
#' @return spec object
calibSpectrum <- function(spec, zoneref, ppmref)
{
   i1<-length(which(spec$ppm>max(zoneref)))
   i2<-which(spec$ppm<=min(zoneref))[1]

   i0 <- i1 + which(spec$int[i1:i2]==max(spec$int[i1:i2])) - 1
   ppm0 <- spec$pmax - (i0-1)*spec$dppm
   dppmref <- ppm0 - ppmref
   decal <- 0
   sig <- C_estime_sd(spec$int[1:spec$size],128)
   if (abs(dppmref) > spec$dppm) {
       decal <- trunc(dppmref/spec$dppm)
       dppmref <- dppmref - decal*spec$dppm
   }
   if (abs(dppmref) > 0.5*spec$dppm) {
       decal <- decal + trunc(2*dppmref/spec$dppm)
       dppmref <- dppmref - trunc(2*dppmref/spec$dppm)*spec$dppm
   }
   if (decal==0) next

   if (decal<0) {
      spec$int[1:(spec$size-abs(decal))] <- spec$int[(1+abs(decal)):spec$size]
      spec$int[(spec$size-abs(decal)+1):spec$size] <- stats::rnorm(length((spec$size-abs(decal)+1):spec$size), mean=spec$int[spec$size-abs(decal)-1], sd=sig)
   }
   if (decal>0) {
      spec$int[(1+abs(decal)):spec$size] <- spec$int[1:(spec$size-abs(decal))]
      spec$int[1:abs(decal)] <- stats::rnorm(length(1:abs(decal)), mean=spec$int[abs(decal)+1], sd=sig)
   }
   return(spec)
}

#' estimateBL
#'
#' \code{estimateBL} estimates of the baseline of the spectrum in the corresponding ppm range (based on the C_Estime_LB routine)
#' @param spec a 'spec' object
#' @param ppmrange the ppm range in which the baseline will be estimated
#' @param WS Size of the window (in number of points) from which a rolling average will be established
#' @param NEIGH The minimum window size (in number of points) in which the signal compared to its mean can be considered as belonging to the baseline. 
#' @return a vector of the estimated baseline
estimateBL <- function(spec, ppmrange, WS=50, NEIGH=35)
{
   set.seed(1234)
   iseq <- c(which(spec$ppm>=ppmrange[1])[1]:length(which(spec$ppm<=ppmrange[2])))
   LB0 <- rep(0,length(iseq))
   # WS : 10, 25, 50
   # NEIGH : 5, 15, 35
   LB0 <- C_Estime_LB (spec$int, min(iseq), max(iseq), WS, NEIGH, 3*spec$Noise)
   LB0
}

#=====================================================================
# GSD - Global Spectra Deconvolution
#=====================================================================

#' sliceSpectrum
#'
#' Slice the spectrum in order to define ranges for Local Deconvolution (LSDeconv)
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param flatwidth specifies the minimum width of a zone in which the spectrum intensities are close to zero to consider this one as a cutting zone (default 0.003 ppm)
#' @param snrfactor specifies factor applied on the Std Dev. of the noise used as threshold for first derivate intensities (default=4)
#' @param maxrange specifies the maximum width of a cutting zone (default 0.3 ppm)
#' @param excludezones specifies the exclusion zones as a matrix (Nx2), each row specifying a zone with 2 columns (ppm min and ppm max) (default NULL)
#' @return a list of ppm range
sliceSpectrum <- function(spec, ppmrange=c(0.5,9.5), flatwidth=0.004, snrfactor=4, maxrange=0.3, excludezones=NULL)
{
   # cut the ppmrange based on first derivate
   # sfac  : factor applied on the Std Dev. of the noise used as threshold for first derivate intensities
   # vfac1 : factor applied on the Std Dev. of the noise used as threshold for spectrum intensities
   # vfac2 : factor applied on the Std Dev. of the noise used as threshold for spectrum intensities
   # delta : min ppm width for flat zones
   cutspectrum <- function(spec, ppmrange, sfac=3, vfac1=20, vfac2=50, delta=0.004) {
      # Parameters
      V <- Rnmr1D:::Smooth(spec$int,10)
      D <- Rnmr1D:::C_Derive1(V)
      SD <- sfac*Rnmr1D:::C_estime_sd(D,64)
      SV <- vfac1*SD
      Vthreshold <- vfac2*SD
      NSIZE <- round(delta/spec$dppm)
      VCUT <- rep(500*SD, length(spec$int))
   
      # Find flat zones
      flg <- nc <- 0
      n1 <- round((ppmrange[1]-spec$pmin)/spec$dppm)
      n2 <- round((ppmrange[2]-spec$pmin)/spec$dppm)
      for (k in n1:n2) {
        if (k==n2) {
          break
        }
        # init
        if (nc==0) {
          nc <- k
          if (abs(D[nc])>SD) flg <- 1
          next
        }
        # No change
        if ( (flg==0 && abs(D[k])<=SD && abs(V[k])<=SV) || (flg==1 && abs(D[k])>SD) ) next
        # Start of a flat zone
        if (flg==1 && abs(D[k])<=SD) {
          nc <- k
          flg <- 0
          next
        }
        # End of the flat zone
        if (flg==0 && (abs(D[k])>SD || abs(V[k])>SV)) {
          if ((k-nc)>=NSIZE) VCUT[nc:(k-1)] <- 0
          nc <- k
          flg <- 1
          next
        }
      }
   
      # Remove false peak zones
      flg <- nc <- 0
      for (k in n1:n2) {
        if ((VCUT[k]>0 & flg==1) || (VCUT[k]==0 & flg==0))
          next
        if (VCUT[k]>0 & flg==0) {
          nc <- k; flg <- 1;
          next
        }
        if ((k-1-nc)<NSIZE || max(V[nc:k-1])<Vthreshold) VCUT[nc:k] <- 0
        flg <- 0
      }
   
      # Get the ppm range of the peaGet the ppm range of the peak zones 
      # by withdrawing the exclude zonesk zones minus the exclude zones
      DN <- round(0.008/spec$dppm)
      flg <- 0
      p0 <- p1 <- NULL
      rangelist <- NULL
      for (k in n1:n2) {
        # No change
        if ((VCUT[k]>0 & flg==1) || (VCUT[k]==0 & flg==0)) next
        # Down to Up
        if (VCUT[k]>0 & flg==0) {
          nup <- k; flg <- 1;
          ncut <- nup - DN
          if (is.null(p0)) {
            p0 <- max(spec$pmin, spec$ppm[ncut])
          } else {
            pcut <- spec$ppm[ncut+which(V[ncut:nup]==min(V[ncut:nup]))-1]
            p0 <- max(spec$pmin, pcut)
          }
          next
        }
        # Up to Down: (VCUT[k]==0 & flg==1)
        ndwn <- k;  flg <- 0
        ncut <- ndwn + DN
        pcut <- spec$ppm[ndwn+which(V[ndwn:ncut]==min(V[ndwn:ncut])) - 1]
        p1 <- min(spec$pmax, pcut)
        fok <- 1
        if (!is.null(excludezones)) {
           for (i in 1:nrow(excludezones)) {
             EZ1 <- excludezones[i,1]; EZ2 <- excludezones[i,2]
             if (p1<EZ1 || p0>EZ2) next
             if (p0<EZ1 && p1>EZ1) {
               if (round((EZ1-p0)/spec$dppm)>NSIZE) {
                 iseq <- getseq(spec,c(p0,EZ1))
                 nm <- which(V[iseq]==max(V[iseq]))
                 if (nm>0.33*length(iseq) && nm <0.66*length(iseq)) rangelist <- rbind(rangelist, c(p0,EZ1))
               }
               p0 <- EZ1
             }
             if (p0<EZ2 && p1>EZ2) {
               if (round((p1-EZ2)/spec$dppm)>NSIZE) {
                 iseq <- getseq(spec,c(EZ2,p1))
                 nm <- which(V[iseq]==max(V[iseq]))
                 if (nm>0.33*length(iseq) && nm <0.66*length(iseq)) rangelist <- rbind(rangelist, c(EZ2,p1))
               }
               p1 <- EZ2
             }
             if (p0>=EZ1 && p1<=EZ2) { fok <- 0; break }
           }
        }
        if (fok)
          rangelist <- rbind(rangelist, c(p0,p1))
      }
      rangelist
   }

   # cut the ppmrange based on the minimum intensities of the spectrum in the ppmrange
   # rangelist : list of ppm ranges to be completed (may be NULL)
   # ppm : ppm range to cut
   # max_ppm_width : max size of ppm ranges as output
   # alpha : proportion of the ppm range to be removed at both ends of the range before searching for the minimums.(0<=alpha<0.5)
   cutrange <- function(spec, rangelist, ppm, max_ppm_width=0.3, alpha=0.25) {
      if (alpha>0.5) alpha <- 1- alpha
      ppmrange <- c( (1-alpha)*ppm[1]+alpha*ppm[2], (1-alpha)*ppm[2]+alpha*ppm[1] )
      iseq <- getseq(spec,ppmrange)
      n <- which(spec$int[iseq]==min(spec$int[iseq]))
      pcut <- spec$ppm[iseq[1]+n-1]
      if ((pcut-ppm[1])<max_ppm_width) {
        rangelist <- rbind( rangelist, c(ppm[1], pcut) )
      }
      else {
        rangelist <- cutrange(spec, rangelist, c(ppm[1], pcut), max_ppm_width)
      }
      if ((ppm[2]-pcut)<max_ppm_width) {
        rangelist <- rbind( rangelist, c(pcut, ppm[2]) )
      }
      else {
        rangelist <- cutrange(spec, rangelist, c(pcut, ppm[2]), max_ppm_width)
      }
      rangelist
   }

   # First cutting 
   rangelist <- cutspectrum(spec, ppmrange, sfac=snrfactor, delta=flatwidth)

   # Second cutting for ppm range greater the specified max ppm range
   rangelist2 <- NULL
   if ("numeric" %in% class(rangelist)) rangelist <- t(as.matrix(rangelist))
   for (k in 1:nrow(rangelist)) {
     ppm <- rangelist[k,]
     if ((ppm[2]-ppm[1])>maxrange) {
       rangelist2 <- cutrange(spec, rangelist2, ppm, maxrange)
     } else {
       rangelist2 <- rbind(rangelist2, ppm)
     }
   }
   rownames(rangelist2) <- NULL
   rangelist2
}

#' getSlices
#'
#' Slice the spectrum in order to define ranges for Local Deconvolution (LSDeconv) and return only those include the provided ppmranges
#' @param spec a 'spec' object
#' @param ppmranges  ppm ranges as a matrix in order to apply the deconvolution, each row specifying a zone
#' @param flatwidth specifies the minimum width of a zone in which the spectrum intensities are close to zero to consider this one as a cutting zone (default 0.003 ppm)
#' @param snrfactor specifies factor applied on the Std Dev. of the noise used as threshold for first derivate intensities (default=4)
#' @param maxrange specifies the maximum width of a cutting zone (default 0.3 ppm)
#' @return a list of ppm range
getSlices <- function(spec, ppmranges, flatwidth=0.004, snrfactor=4, maxrange=0.3)
{
  # Parameters
  PPMRANGE <- 0.9*c(spec$pmin, spec$pmax)
  
  # get the ppm slice list
  rangelist <- Rnmr1D::sliceSpectrum(spec, PPMRANGE, flatwidth=flatwidth, snrfactor=snrfactor, maxrange=maxrange)
  if ("numeric" %in% class(rangelist)) rangelist <- t(as.matrix(rangelist))

  # Range list selection
  V <- NULL
  if ("numeric" %in% class(ppmranges)) ppmranges <- t(as.matrix(ppmranges))
  V <- NULL
  for (p in 1:nrow(ppmranges)) {
    for (r in 1:nrow(rangelist)) {
      P1 <- ppmranges[p,1]; P2 <- ppmranges[p,2]; 
      R1 <- rangelist[r,1]; R2 <- rangelist[r,2];
      if (R1>P2 || R2<P1) next
      R1 <- min(R1,P1)
      R2 <- max(R2,P2)
      V <- rbind(V, c(R1,R2))
      break
    }
  }
  # Overlapping
  for (k in 2:nrow(V)) {
    if (V[k,1]>V[k-1,2]) next
    n1 <- round((V[k,1]-spec$pmin)/spec$dppm)
    n2 <- round((V[k-1,2]-spec$pmin)/spec$dppm)
    vm <- ifelse( spec$int[n1] < spec$int[n2], V[k,1], V[k-1,2] )
    V[k,1] <- V[k-1,2] <- vm
  }
  unique(V[order(V[,1]), ])
}

getBestPeaks <- function(spec, model1, model2, crit=0)
{
   Peaks <- NULL
   nblk1 <- model1$blocks$cnt
   nblk2 <- model2$blocks$cnt
   for (k in 1:nblk1) {
      n1 <- model1$blocks$blocks[k,1]
      n2 <- model1$blocks$blocks[k,2]
      Y0 <- spec$int[n1:n2]
      Y1 <- model1$model[n1:n2]; r1 <- stats::cor(Y0,Y1)^2; sd1 <- stats::sd(Y0-Y1)
      Y2 <- model2$model[n1:n2]; r2 <- stats::cor(Y0,Y2)^2; sd2 <- stats::sd(Y0-Y2)
      if (crit==0) { c <- (r1>r2) } else { c <- (sd1<sd2) }
      if (c) { P2 <- model1$peaks[model1$peaks$pos>n1, ] }
      else   { P2 <- model2$peaks[model2$peaks$pos>n1, ] }
      P2 <- P2[ P2$pos<n2, ]
      Peaks <- rbind(Peaks, P2)
   }
   Peaks
}

#' GSDeconv
#'
#' Global Spectra Deconvolution: \code{GSDeconv} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param filter a filter type for filtering the noise and  smoothing the signal
#' @param scset a set of scmin values
#' @param params a list of specific parameters for deconvolution
#' @param verbose level of debug information
#' @return a model object
GSDeconv <- function(spec, ppmrange, params=NULL, filter='symlet8', scset=c(2,3,12), verbose=1)
{

   if ( ! 'Spectrum' %in% class(spec) )
      stop("the input spec must belong to the 'Spectrum' class")
   g <- getDeconvParams(params)
   g$oneblk <- 0;

   FacN <- ifelse( is.null(g$facN), 5, g$facN )
   spec$B <- spec$Noise/FacN
   g$ratioSN <- ifelse(g$lowPeaks==0, FacN*g$ratioPN, FacN/10)
   debug1 <- ifelse(verbose==2, 1, 0)

   # Peak search
   model0 <- C_peakFinder(spec, ppmrange, g$flist[[filter]], g, verbose=0)

   # Peak optimization
   g$peaks <- model0$peaks
   vset <- scset[order(scset)]
   for (k in vset) {
      if (debug1) cat(k,':')
      g$scmin <- k
      if (k==min(vset)) {
         model1 <- C_peakOptimize(spec, ppmrange, g, verbose = debug1)
      } else {
         model2 <- C_peakOptimize(spec, ppmrange, g, verbose = debug1)
         model1$peaks <- getBestPeaks(spec, model1, model2, crit=g$criterion)
      }
   }
   Ymodel <- specModel(spec, ppmrange, model1$peaks)
   iseq <- getIndexSeq(spec,ppmrange)

   model1$model <- Ymodel
   model1$R2 <- stats::cor(spec$int[iseq],Ymodel[iseq])^2
   model1$residus <- spec$int-Ymodel
   model1$SD <- stats::sd(model1$residus[iseq]/spec$Noise)

   if (verbose) {
      cat('FacN =',FacN,', RatioPN =',g$ratioPN,', RatioSN =',g$ratioPN*FacN,"\n")
      cat('Nb Blocks =',model1$blocks$cnt,',Nb Peaks =', model1$nbpeak,"\n")
      cat('R2 =', model1$R2,"\n")
      cat('Residue/Noise : SD =',round(model1$SD,4), ', Mean =',round(mean(model1$residus)/spec$Noise,4), "\n")
   }
   class(model1) = "GSDmodel"
   model1
}

#=====================================================================
# LSD - Local Spectra Deconvolution
#=====================================================================

intern_computeBL <- function(spec, model)
{
  repeat {
     bl <- rep(0,length(spec$ppm))
     if (model$params$obl==0) break

     nblk <- model$blocks$cnt
     if (nblk==0) break

     blocks <- model$blocks$blocks
     blpars <- model$blocks$blpars

     for (k in 1:nblk) {
        iseq <- blocks[k,1]:blocks[k,2]
        pm <- (blocks[k,3]+blocks[k,4])/2
        a <- blpars[k, ]
        p <- spec$ppm[iseq]
        xp <- 1;
        for (k in 1:length(a)) { bl[iseq] <- bl[iseq] + a[k]*xp; xp <- xp*(p - pm); }
     }
     break
  }
  bl
}

#' computeBL
#'
#' \code{computeBL} computes baseline based on the model.
#' @param spec a 'spec' object
#' @param model a model object
#' @return a vector of the baseline estimated during the deconvolution process
computeBL <- function(spec, model)
{
   if ( ! 'Spectrum' %in% class(spec) )
      stop("the input spec must belong to the 'Spectrum' class")
   if ( ! sum(c('optimModel', 'LSDmodel') %in% class(model) ) )
      stop("the input model must have an appropriate class, namely 'optimModel' or 'LSDmodel'")
   intern_computeBL(spec, model)
}

#' LSDeconv
#'
#' Local Spectra Deconvolution: \code{LSDeconv} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param params a list of specific parameters for deconvolution including or not (i.e equal to NULL) the matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, eta
#' @param filterset a set of filter type for filtering the noise and  smoothing the signal (only if the matrix defining peaks not defined in order to find peaks)
#' @param oblset a set of baseline order for fitting
#' @param verbose level of debug information
#' @return a model object
LSDeconv <- function(spec, ppmrange, params=NULL, filterset=1:6, oblset=1:12, verbose=1)
{
   if ( ! 'Spectrum' %in% class(spec) )
      stop("the input spec must belong to the 'Spectrum' class")

   set.seed(1234)
   if (is.null(params$peaks)) {
      model <- LSDeconv_1(spec, ppmrange, params, filterset, oblset, verbose)
      if (params$spass>0) { # second deconvolution phase without the highest peaks has to be done
          model <- LSDspass(spec, model, ppmrange, params, oblset, verbose)
      }
   } else {
      model <- LSDeconv_2(spec, ppmrange, params, oblset, verbose)
   }
   model
}

#' LSDspass
#'
#' Local Spectra Second Deconvolution Phase: \code{LSDspass} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param model a 'model' object 
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param params a list of specific parameters for deconvolution including or not (i.e equal to NULL) the matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, eta
#' @param oblset a set of baseline order for fitting
#' @param verbose level of debug information
#' @return a model object
LSDspass <- function(spec, model, ppmrange, params=NULL, oblset=1:12, verbose=1)
{
   # Select significant peaks
   pk1 <- model$peaks[model$peaks$amp/spec$Noise>100, ]
   # if some significant peaks have an intensity less than 20% of the highest peak
   if (nrow(pk1)>0 && (pk1$amp/max(pk1$amp))<0.2) {
       # Select significant peaks with intensity greater than 50% of the highest peak
       hpkpos <- pk1[(pk1$amp/max(pk1$amp))>0.5, ]$pos
       if (length(hpkpos)<model$nbpeak) {
          pksel <- model$peaks[which(model$peaks$pos == hpkpos), ]
          Y1 <- specModel(spec, ppmrange, pksel)
          specInt <- spec$int
          spec$int <- specInt - Y1
          model_filter <- model$filter
          # Select the other peaks to be deconvoluted separately
          params$peaks <- model$peaks[which(model$peaks$pos != hpkpos), ]
          model <- LSDeconv_2(spec, ppmrange, params, oblset, verbose = verbose)
          # Add the previous selected peaks in the model
          spec$int <- specInt
          model$filter <- model_filter
          model$peaks <- rbind(model$peaks, pksel)
          model$peaks <- model$peaks[order(model$peaks$pos),]
          model$nbpeak <- model$nbpeak + nrow(pksel)
          rownames(model$peaks) <- 1:model$nbpeak
          model$model <- specModel(spec, ppmrange, model$peaks)
          model$residus <- spec$int - model$model - model$LB
          model$RMSE <- sqrt(mean(model$residus^2))
          model$R2 <- stats::cor(spec$int[model$iseq],model$model[model$iseq])^2
       }
   }
   model
}

# Local Spectra Deconvolution with no predefined peaks (g$peaks=NULL)
LSDeconv_1 <- function(spec, ppmrange, params=NULL, filterset=1:6, oblset=1:12, verbose=1)
{
   g <- getDeconvParams(params)
   #g$oneblk <- 1;
   iseq <- getIndexSeq(spec,ppmrange)
   if (is.null(g$facN))
      FacN <- max(20,min(100,round(max(spec$int[iseq])*0.05/spec$Noise)))
   else
      FacN <- g$facN
   spec$B <- spec$Noise/FacN
   g$ratioSN <- ifelse(g$lowPeaks==0, FacN*g$ratioPN, FacN/10)

   OBL <- NULL
   R2 <- NULL
   SD <- NULL
   debug1 <- ifelse(verbose==2, 1, 0)

# Set of values for filter
   for (filt in filterset) {
      R2i <- NULL
      SDi <- NULL
# Set of values for order of the polynomial model of the baseline
      for (obl in oblset) {
         g$obl <- obl
         model0 <- C_peakFinder(spec, ppmrange, g$flist[[filt]], g, verbose = 0)
         model0$peaks <- Rnmr1D::peakFiltering(spec,model0$peaks, g$ratioSN)
         model0$nbpeak <- ifelse('data.frame' %in% class(model0$peaks), nrow(model0$peaks), 0)
         if (model0$nbpeak==0) {
            R2i <- c( R2i, 0 ); SDi <- c( SDi, 1e999 );
            next
         }
         g$peaks <- model0$peaks
         model <- C_peakOptimize(spec, ppmrange, g, verbose = 0)
         model$peaks <- Rnmr1D::peakFiltering(spec,model$peaks, g$ratioPN*FacN)
         model$model <- Rnmr1D::specModel(spec, ppmrange, model$peaks)
         if (model$nbpeak>0) {
            Ymodel <- model$model + intern_computeBL(spec, model)
            R2i <- c( R2i, stats::cor(spec$int[iseq],Ymodel[iseq])^2 )
            SDi <- c( SDi, model$blocks$blocks[1,6] )
         } else {
            R2i <- c( R2i, 0 ); SDi <- c( SDi, 1e999 )
         }
         idx <- length(R2i)
if (debug1) cat(filt,": R2i =",round(R2i[idx],4),", RMSEi =",round(SDi[idx],6),", PKi =",model0$nbpeak," - ",round(mean(model0$peaks$sigma),6),
                     " OBL =",obl," ETA =",model$peaks$eta[1],"\n")
      }
      idx <- ifelse ( g$criterion==0, which(R2i==max(R2i))[1], which(SDi==min(SDi))[1] )
      OBL <- c( OBL, oblset[idx] )
      R2 <- c( R2, R2i[idx] )
      SD <- c( SD, SDi[idx] )
if (debug1) cat(filt,": R2 =",round(R2i[idx],4),", RMSE =",round(SDi[idx],6)," OBL =",oblset[idx],"\n\n")
   }
   gc()

   R2[ is.na(R2) ] <- 0; SD[ is.na(SD) ] <- 1e999
   if (is.null(R2) || sum(R2)==0)
      stop("No peak found.")

   idx <- ifelse ( g$criterion==0, which(R2==max(R2)), which(SD==min(SD)) )
   fidx <- filterset[idx]
   if (debug1) cat("Best: idx =",idx,", filter =",fidx,", obl =",OBL[idx],"\n")

   model0 <- C_peakFinder(spec, ppmrange, g$flist[[fidx]], g, verbose = debug1)
   if (debug1) cat("----\n")
   g$obl <- OBL[idx]
   g$peaks <- model0$peaks
   model <- C_peakOptimize(spec, ppmrange, g, verbose = debug1)
   model$peaks <- Rnmr1D::peakFiltering(spec,model$peaks, g$ratioPN*FacN)
   model$model <- Rnmr1D::specModel(spec, ppmrange, model$peaks)

   if (debug1) cat("----\n")
   P1 <- model$peaks[model$peaks$amp>0, ]
   P2 <- P1[P1$ppm>ppmrange[1], ]
   model$peaks <- P2[P2$ppm<ppmrange[2],]
   rownames(model$peaks) <- NULL
   model$nbpeak <- nrow(model$peaks)
   model$LB <- intern_computeBL(spec, model)
   Ymodel <- model$model + model$LB

   model$residus <- spec$int-Ymodel
   model$iseq <- iseq
   model$ppmrange <- ppmrange
   model$R2 <- stats::cor(spec$int[iseq],Ymodel[iseq])^2
   model$filter <- fidx
   model$crit <- g$crit
   model$FacN <- FacN
   model$eta <- mean(model$peaks$eta)
   model$RMSE <- sqrt(mean((model$residus[iseq])^2))

   if (verbose) {
      cat('FacN =',FacN,', RatioPN =',g$ratioPN,', RatioSN =',g$ratioSN,"\n")
      cat('crit =',model$crit,', filter =', model$filter,', obl =',model$params$obl,', eta =',model$eta,"\n")
      cat('Nb Blocks =',model$blocks$cnt,', Nb Peaks =', model$nbpeak,"\n")
      cat('RMSE =', model$RMSE,"\n")
      cat('R2 =', model$R2,"\n")
      cat('Residue : SD =',round(stats::sd(model$residus[iseq]),4),
                ', Mean =',round(mean(model$residus[iseq]),4), "\n")
   }
   class(model) = "LSDmodel"
   model
}

# Local Spectra Deconvolution with predefined peaks (g$peaks not NULL)
LSDeconv_2 <- function(spec, ppmrange, params=NULL, oblset=1:12, verbose=1)
{
   g <- getDeconvParams(params)
   #g$oneblk <- 1;
   iseq <- getIndexSeq(spec,ppmrange)
   if (is.null(g$facN))
      FacN <- max(20,min(100,round(max(spec$int[iseq])*0.05/spec$Noise)))
   else
      FacN <- g$facN
   spec$B <- spec$Noise/FacN
   g$ratioSN <- ifelse(g$lowPeaks==0, FacN*g$ratioPN, FacN/10)

   if (is.null(g$peaks) || ! "data.frame" %in% class(g$peaks) || nrow(g$peaks)==0 )
      stop("the peaks param must be a data.frame with at least one row")

   R2 <- NULL
   SD <- NULL
   debug1 <- ifelse(verbose==2, 1, 0)
   for (obl in oblset) {
      g$obl <- obl
      model <- C_peakOptimize(spec, ppmrange, g, verbose = debug1)
      if (model$nbpeak>0) {
         Ymodel <- model$model + intern_computeBL(spec, model)
         residus <- spec$int-Ymodel
         R2 <- c( R2, stats::cor(spec$int[iseq],Ymodel[iseq])^2 )
         SD <- c( SD, stats::sd(residus[iseq]/spec$Noise) )
      } else {
         R2 <- c( R2, 0 ); SD <- c( SD, 1e999 )
      }
   }
   gc()

   R2[ is.na(R2) ] <- 0; SD[ is.na(SD) ] <- 1e999
   if (is.null(R2) || sum(R2)==0)
      stop("No peak found.")

   idx <- ifelse ( g$criterion==0, which(R2==max(R2))[1], which(SD==min(SD))[1] )
   if (debug1) cat("Best: idx =",idx,", obl =",oblset[idx],"\n")

   g$obl <- oblset[idx];
   model <- C_peakOptimize(spec, ppmrange, g, verbose = debug1)
   P1 <- model$peaks[model$peaks$ppm>ppmrange[1], ]
   model$peaks <- P1[P1$ppm<ppmrange[2],]
   rownames(model$peaks) <- NULL
   model$nbpeak <- dim(model$peaks)[1]
   model$LB <- intern_computeBL(spec, model)
   Ymodel <- model$model + model$LB

   model$residus <- spec$int-Ymodel
   model$iseq <- iseq
   model$ppmrange <- ppmrange
   model$R2 <- stats::cor(spec$int[iseq],Ymodel[iseq])^2
   model$crit <- g$crit
   model$FacN <- FacN
   model$RMSE <- sqrt(mean(model$residus^2))

   if (verbose) {
      cat('FacN =',FacN,', RatioPN =',g$ratioPN,', RatioSN =',g$ratioPN*FacN,"\n")
      cat('crit =',model$crit,', obl =',model$params$obl,', eta =',model$peaks$eta[1],"\n")
      cat('Nb Blocks =',model$blocks$cnt,', Nb Peaks =', model$nbpeak,"\n")
      cat('RMSE =', model$RMSE,"\n")
      cat('R2 =', model$R2,"\n")
      cat('Residue : SD =',round(stats::sd(model$residus[iseq]),4),
                ', Mean =',round(mean(model$residus[iseq]),4), "\n")
   }
   class(model) = "LSDmodel"
   model
}

#' MultiLSDeconv
#'
#' Multiple Local Spectra Deconvolution: \code{MultiLSDeconv} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmranges ppm ranges as a matrix in order to apply the deconvolution, each row specifying a zone
#' @param params a list of specific parameters for deconvolution 
#' @param filterset a set of filter type for filtering the noise and  smoothing the signal (only if the matrix defining peaks not defined in order to find peaks)
#' @param oblset a set of baseline order for fitting
#' @param ncpu number of CPU for parallel computing
#' @param verbose level of debug information
#' @return a model object
MultiLSDeconv <- function(spec, ppmranges=NULL, params=NULL, filter='symlet8', oblset=0, ncpu=4, verbose=0)
{
   g <- getDeconvParams(params)
   excludezones <- g$exclude_zones
   if (is.null(ppmranges)) {
       slices <- sliceSpectrum(spec, flatwidth=g$flatwidth, snrfactor=g$snrfactor, maxrange=g$maxrange, excludezones=excludezones)
   } else {
       slices <- getSlices(spec, ppmranges, flatwidth=g$flatwidth, snrfactor=g$snrfactor, maxrange=g$maxrange)
   }
   if (verbose) cat("Nb slices =",nrow(slices),", minsize =",min(slices[,2]-slices[,1]),", maxsize =",max(slices[,2]-slices[,1]),"\n")

   combine_list <- function(LL1, LL2) {
      getList <- function(L) { list(id=L$id, peaks=L$peaks, infos=L$infos, LB=L$LB ) }
      mylist <- list()
      for ( i in 1:length(LL1) ) mylist[[LL1[[i]]$id]] <- getList(LL1[[i]])
      for ( i in 1:length(LL2) ) mylist[[LL2[[i]]$id]] <- getList(LL2[[i]])
      return(mylist)
   }

   # Activate the cluster
   cl <- parallel::makeCluster(ncpu)
   doParallel::registerDoParallel(cl)

   M <- foreach( k = 1:nrow(slices), .combine = combine_list, .packages=c('Rnmr1D'), .export=c('getIndexSeq') ) %dopar% {
      iseq <- getIndexSeq(spec,slices[k,])
      facN <- max(20,min(100,round(max(spec$int[iseq])*0.05/spec$Noise)))
      spec$B <- spec$Noise/facN
      g$ratioSN <- facN*g$ratioPN
      debug1 <- ifelse(verbose>1, 1, 0)
      t<-system.time({
         model <- tryCatch({
             LSDeconv(spec, slices[k,], g, filter, oblset, verbose = debug1)
         },
         error=function(e) {
             model <- list(peaks=NULL, LB=NULL)
         })
      })
      infos <- NULL
      if (!is.null(model$peaks)) {
         infos <- c(round(slices[k,1],4), round(slices[k,2],4),round(slices[k,2]-slices[k,1],4), model$nbpeak, 
                    round(model$eta,4), model$filter, model$params$obl, round(model$R2,5), round(t[3],2))
      }
      id <- paste0('S',sprintf("%03d",k))
      mylist <- list()
      mylist[[id]] <- list(id=id, peaks=model$peaks, infos=infos, LB=model$LB)
      return(mylist)
   }

   # Stop the cluster
   parallel::stopCluster(cl)

   peaks <- infos <- NULL
   LB <- rep(0,length(spec$ppm))
   for (k in 1:nrow(slices)) {
      id <- paste0('S',sprintf("%03d",k))
      if (!is.null(M[[id]]$peaks)) {
         peaks <- rbind(peaks, M[[id]]$peaks)
         infos <- rbind(infos, M[[id]]$infos)
         iseq <- getIndexSeq(spec,slices[k,])
         LB[iseq] <- M[[id]]$LB[iseq]
      }
   }
   colnames(infos) <- c("ppm1","ppm2","width","peaks","eta","filter","obl","R2","time")
   rownames(infos) <- 1:nrow(infos)

   ppmrange <- c(spec$pmin, spec$pmax)
   Ymodel <- specModel(spec, ppmrange, peaks)
   if (is.null(ppmranges)) {
      yorig <- rep(0, length(spec$int))
      iseq <- getIndexSeq(spec,ppmrange)
      yorig[iseq] <- spec$int[iseq]
      if (! is.null(excludezones) && nrow(excludezones)>0) {
         for (z in 1:nrow(excludezones))
            yorig[getIndexSeq(spec,excludezones[z,])] <- 0
      }
   } else {
      yorig <- rep(0, length(spec$int))
      for (k in 1:nrow(slices)) {
          if (!is.null(M[[id]]$peaks)) {
             iseq <- getIndexSeq(spec,slices[k,])
             yorig[iseq] <- spec$int[iseq]
          }
      }
   }
   residus <- yorig - Ymodel - LB
   RMSE <- sqrt(mean(residus^2))
   R2 <- stats::cor(yorig,Ymodel+LB)^2
   model <- list(peaks=peaks, infos=infos, slices=slices, model=Ymodel, LB=LB, residus=residus, R2=R2, RMSE=RMSE, filter=filter, params=g)
   model
}

#' cleanPeaks
#'
#' \code{cleanPeaks} cleans the peaks under a specified threshold 
#' @param spec a 'spec' object
#' @param model a model object
#' @param SNthreshold Threshold for the Signal-Noise Ratio below which the peaks will be rejected
#' @return a data.frame of the remaining peaks
cleanPeaks <- function(spec, model, SNthreshold=5) {
   cmodel <- model$peaks[model$peaks$amp/spec$Noise > SNthreshold, -5 ]
   ppm <- model$ppmrange
   cmodel <- cmodel[cmodel$ppm>ppm[1],]; cmodel <- cmodel[cmodel$ppm<ppm[2],]
   cmodel$PNratio<- cmodel$amp/spec$Noise
   rownames(cmodel) <- c(1:length(cmodel$ppm))
   cmodel
}

#=====================================================================
# Plots
#=====================================================================


#' plotSpec
#'
#' \code{plotSpec} plots all signals defined by a matrix.
#' @param ppmrange a ppm range defining the window of the plotting
#' @param x a vector defining the x-axis (abscissa)
#' @param y a vector or a matrix defining the y-axes (ordinates), each signal as a column
#' @param ynames a vector defining the y names (same order as the y matrix)
#' @param ycolors a vector defining the y colors (same order as the y matrix)
#' @param ysel a vector defining the visibility of each y element (same order as the y matrix)
#' @param title title of the graphic
plotSpec <- function(ppmrange, x, y, ynames=c('Origin', 'Filtered', 'Model'), 
                     ycolors=c('grey', 'blue', 'red', 'green', 'orange','magenta','cyan','darkgreen', 'darkorange'), ysel=NULL, title='')
{
   iseq <- c(which(x>=ppmrange[1])[1]:length(which(x<=ppmrange[2])))
   if ("numeric" %in% class(y)) {
      data <- data.frame(x=x[iseq], y=y[iseq])
      p <- plotly::plot_ly(data, x = ~x, y = ~y, name = ynames[1], type = 'scatter', mode = 'lines')
   }
   else {
      if (is.null(ysel)) ysel <- rep(TRUE,ncol(y))
      y1 <- y[,1]
      data <- data.frame(x=x[iseq], y=y1[iseq])
      visible <- ifelse(ysel[1],  TRUE , "legendonly" )
      p <- plotly::plot_ly(data, x = ~x, y = ~y, name = ynames[1], type = 'scatter', mode = 'lines', visible=visible)
      for (k in 2:ncol(y)) {
         df <- data.frame(x=x[iseq], y=y[iseq,k])
         visible <- ifelse(ysel[k],  TRUE , "legendonly" )
         p <- p %>% plotly::add_trace(data=df, x = ~x, y = ~y, name=ynames[k], mode = 'lines', visible=visible)
      }
   }
   p <- p %>% plotly::layout(title=title, xaxis = list(autorange = "reversed"), colorway = ycolors)
   p
}

#' plotModelwithResidus
#'
#' \code{plotSpec} plots the model along with the resulting residues from deconvolution
#' @param spec a 'spec' object (see \code{readSpectrum}, \code{Spec1rDoProc})
#' @param model a 'model' object (see \code{specDeconv}, \code{peakOptimize}, \code{GSDeconv}, \code{LSDeconv})
#' @param ynames a vector defining the y names. default value = c('Origin', 'Model', 'Residues')
#' @param title title of the graphic
plotModelwithResidus <- function(spec, model, ynames=c('Origin', 'Model', 'Residues'), title='')
{
   if ( ! sum(c('GSDmodel', 'LSDmodel') %in% class(model) ) )
      stop("the input model must have an appropriate class, namely 'GSDmodel' or 'LSDmodel'")
   Ymodel <- model$model
   if (model$params$obl>0) Ymodel <- Ymodel + model$LB
   iseq <- model$iseq
   data <- data.frame( x=spec$ppm[iseq], y=spec$int[iseq], y1=Ymodel[iseq], y2=model$residus[iseq] )
   p <- plotly::plot_ly(data, x = ~x, y = ~y, name = ynames[1], type = 'scatter', mode = 'lines', 
                fill = 'tozeroy', fillcolor='rgba(220,220,220,0.1)') %>%
     plotly::add_trace(y = ~y1, name = ynames[2], mode = 'lines', hoverinfo = 'skip') %>%
     plotly::add_trace(y = ~y2, name = ynames[3], mode = 'lines') %>%
     plotly::layout(title=title, xaxis = list(autorange = "reversed"), colorway = c('grey', 'blue', 'red'))
   p
}

#' plotModel
#'
#' \code{plotModel} plots the model along with the resulting voigt functions from deconvolution
#' @param spec a 'spec' object (see \code{readSpectrum}, \code{Spec1rDoProc})
#' @param model a 'model' object (see \code{specDeconv}, \code{peakOptimize}, \code{GSDeconv}, \code{LSDeconv})
#' @param exclude_zones a list of vector defining the excluded zones for lorentzian plots
#' @param labels choose as legend labels either 'ppm' or 'id'
#' @param tags boolean allowing you to put identification tags at the top of each peak
#' @param title title of the graphic
plotModel <- function(spec, model, exclude_zones=NULL, labels=c('ppm','id'), tags=FALSE, title='')
{
   if ( ! sum(c('peakModel', 'optimModel','GSDmodel', 'LSDmodel') %in% class(model) ) )
      stop("the input model must have an appropriate class, namely one of these: 'peakModel', 'optimModel', 'GSDmodel', 'LSDmodel'")
   ppmrange <- c(model$params$wmin, model$params$wmax)
   iseq <- getIndexSeq(spec,ppmrange)
   P1 <- model$peaks[model$peaks$ppm>ppmrange[1], ]
   P2 <- P1[P1$ppm<ppmrange[2],]
   if (! is.null(exclude_zones))
     for (k in 1:length(exclude_zones)) P2 <- rbind( P2[P2[,2]<exclude_zones[[k]][1],], P2[P2[,2]>exclude_zones[[k]][2],] )
   npk <- nrow(P2)
   npk_colors <- sample(grDevices::rainbow(npk, s=0.8, v=0.75))
   V <- simplify2array(lapply(1:npk, function(i) { PVoigt(spec$ppm[iseq], P2$amp[i], P2$ppm[i], P2$sigma[i], P2$asym[i], P2$eta[i]) }))
   fmodel <- apply(V,1,sum)
   datamodel <- data.frame(x=spec$ppm[iseq], ymodel=fmodel)
   labels <- match.arg(labels)
   p1 <- plotly::plot_ly(datamodel, x = ~x, y = ~ymodel, name = 'Model', type = 'scatter', mode = 'lines' )
   for (i in 1:npk){
     df <- data.frame(x=datamodel$x, y=V[,i])
     p1 <- plotly::add_trace(p1, data=df, x = ~x, y = ~y, name=ifelse(labels=='ppm', paste0("p",round(P2[i,2],5)), i), mode = 'lines', fill = 'tozeroy')
   }
   p1 <- p1 %>% plotly::layout(title = title, xaxis = list(autorange = "reversed"), colorway = c('grey', npk_colors))
   if (tags) {
     data <- data.frame(lab=rownames(P2), x=P2[,2], y=1.05*P2[,3])
     p1 <- p1 %>% plotly::add_annotations(x = data$x, y = data$y, text = as.character(data$lab), showarrow = TRUE, arrowcolor='red', 
                           font = list(color = 'black', family = 'sans serif', size = 18))
   }
   p1
}
