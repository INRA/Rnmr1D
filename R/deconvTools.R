#------------------------------------------------
# Rnmr1D package: ID deconvTools.R
# Project: NMRProcFlow
# (C) 2015-2021 - D. JACOB - IMRAE UMR1332 BAP
#------------------------------------------------

#=====================================================================
# Deconvolution parameters
#=====================================================================


# Filter types
fnone <- list( type=0 )
fdaub8 <- list( type=1, threshold=0.5 )
fsymlet8 <- list( type=2, threshold=0.5 )
fsavgol5 <- list( type=3, m=5, nl=8, nr=8 )
fsavgol10 <- list( type=3, m=10, nl=20, nr=30 )
fsavgol50 <- list( type=3, m=50, nl=60, nr=80 )
fsmooth <- list( type=4, m=12 )

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
#'   \item \code{dist_fac} : PeakFinder : min distance between 2 peaks (as multiple of sigma_min which is typically equal to 0.0005 ppm) - default value = 2
#'   \item \code{optim} : Indicates if optimisation is applied - default value = 1
#'   \item \code{oppm} : Indicates if ppm optimisation is applied - default value = 1
#'   \item \code{osigma} : Indicates if sigma optimisation is applied - default value = 1
#'   \item \code{d2meth} : PeakFinder : Indicates if minima method to the second derivation is applied
#'   \item \code{spcv} : PeakFinder : Maximal CV on Spectrum - default value = 0.005
#'   \item \code{d2cv} : PeakFinder : Maximum CV on the derived spectrum - default value = 0.05
#'   \item \code{d1filt} : Apply Filter (1) on the 1st derivate or not (0) - default value = 0
#'   \item \code{d2filt} : Apply Filter (1) on the 2nd derivate or not (0) - default value = 1
#'   \item \code{sigma_min} : Optimization of Sigmas : Fixe the minimum limit of sigmas - default value = 0.0005
#'   \item \code{sigma_max} : Optimization of Sigmas : Fixe the maximum limit of sigmas - default value = 0.005
#'   \item \code{verbose} : Indicates if we want information messages - default value = 1
#'   \item \code{exclude_zones} : Exclude ppm zones for the criterion evaluation - default value = NULL
#' }
deconvParams <- list (
  # Filter types
  flist = list( 'smooth1'=fsavgol5, 'smooth2'=fsavgol10, 'smooth3'=fsavgol50,
                'daub8'=fdaub8, 'symlet8'=fsymlet8, 'none'=fnone ),

  # Criterion type for the optimizations
  # 0 => R2
  # 1 => 1/Std(residues)
  criterion = 0,

  # Criterion tolerance for the optimization
  reltol = 0.0001,

  # Peak/Noise Ratio
  facN = NULL,
  ratioPN = 5,

  # indicates if pseudo-voigt is used instead of lorentzian
  pvoigt=0,

  # Optimization of peaks : 0 => No, 1 => Yes
  optim=1,
  oppm = 1,
  osigma = 1,
  oeta=0,
  estimate_int=0,

  # Optimization by only one block or by several blocks applying a cut-off process. 
  oneblk=0,
  scmin=2,

  # Optimization of a baseline (BL) for each massif
  # 0 means no BL, an integer greater than 0 indicates the polynomial order of the BL
  obl = 0,

  # Peaks searching : min distance between 2 peaks (as multiple of sigma_min which is typically equal to 0.0005 ppm)
  dist_fac = 2,

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
  
  # a dataframe of peaks (columns : pos, ppm, amp, sigma, pfac)
  peaks = NULL
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
#' @return a vector of the lorentzian values (same size as ppm)
Lorentz <- function(ppm, amp, x0, sigma)
{
   C_Lorentz(ppm, amp, x0, sigma)
}

#' PVoigt
#'
#' \code{PVoigt} belongs to the low-level functions group for deconvolution.
#' @param ppm a vector of ppm values
#' @param amp amplitude of the lorentzian
#' @param x0 central value of the lorentzian
#' @param sigma half-width of the lorentzian
#' @param sigma2 half-width of the gaussian
#' @param eta mixing coefficient for the pseudo-voigt function (between 0 and 1)
#' @return a vector of the lorentzian values (same size as ppm)
PVoigt <- function(ppm, amp, x0, sigma, sigma2, eta)
{
   C_PVoigt(ppm, amp, x0, sigma, sigma2, eta)
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
#' @param params a list of specific parameters for deconvolution, including the matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, pfac
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

#' specModel
#'
#' \code{specModel} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param peaks a matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, pfac, integral
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
   g$ratioSN <- FacN*g$ratioPN
   model0 <- peakFinder(spec, ppmrange, g, 'daub8', verbose = 0)
   model0$nbpeak
}

#=====================================================================
# Processing
#=====================================================================

# Read then process the spectrum

# directory containing the FID

#' read_spectrum
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
   spec <- Spec1rDoProc(Input=ACQDIR,param=procParams)
   Noise <- C_noise_estimation(spec$int, which(spec$ppm>=ppmnoise[1])[1], length(which(spec$ppm<=ppmnoise[2])))
   spec$int <- spec$int/scaleIntensity
   spec$B <- spec$Noise <- Noise/scaleIntensity
   if (verbose) {
      cat('Size =',length(spec$int),', Max Intensity =',round(max(spec$int),3),', Noise Level =',round(spec$Noise,6),"\n")
   }
   spec
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
   g <- getDeconvParams(params)
   g$oneblk <- 0;

   FacN <- ifelse( is.null(g$facN), 5, g$facN )
   spec$B <- spec$Noise/FacN
   g$ratioSN <- FacN*g$ratioPN
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
        a <- blpars[k, ]
        pm <- (blocks[k,3]+blocks[k,4])/2
        xp <- 1;
        if (nblk==1) {
            p <- spec$ppm
            for (k in 1:length(a)) { bl <- bl + a[k]*xp; xp <- xp*(p - pm); }
        } else {
            iseq <- blocks[k,1]:blocks[k,2]
            p <- spec$ppm[iseq]
            for (k in 1:length(a)) { bl[iseq] <- bl[iseq] + a[k]*xp; xp <- xp*(p - pm); }
        }
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
   if ( ! sum(c('optimModel', 'LSDmodel') %in% class(model) ) )
      stop("the input model must have an appropriate class, namely 'optimModel' or 'LSDmodel'")
   intern_computeBL(spec, model)
}

#' LSDeconv
#'
#' Local Spectra Deconvolution: \code{LSDeconv} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param params a list of specific parameters for deconvolution including or not (i.e equal to NULL) the matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, pfac
#' @param filterset a set of filter type for filtering the noise and  smoothing the signal (only if the matrix defining peaks not defined in order to find peaks)
#' @param oblset a set of baseline order for fitting
#' @param filterset a set of filter applied before fitting
#' @param verbose level of debug information
#' @return a model object
LSDeconv <- function(spec, ppmrange, params=NULL, filterset=1:6, oblset=1:12, verbose=1)
{
   set.seed(1234)
   if (is.null(params$peaks))
      LSDeconv_1(spec, ppmrange, params, filterset, oblset, verbose)
   else
      LSDeconv_2(spec, ppmrange, params, oblset, verbose)
}

# Local Spectra Deconvolution with no predefined peaks (g$peaks=NULL)
LSDeconv_1 <- function(spec, ppmrange, params=NULL, filterset=1:6, oblset=1:12, verbose=1)
{
   g <- getDeconvParams(params)
   iseq <- getIndexSeq(spec,ppmrange)
   if (is.null(g$facN))
      FacN <- max(20,min(100,round(max(spec$int[iseq])*0.05/spec$Noise)))
   else
      FacN <- g$facN
   spec$B <- spec$Noise/FacN
   g$ratioSN <- FacN*g$ratioPN

   etaset <- c(g$eta)
   if (g$oeta==1) etaset <-seq(0.60,0.75,0.025)
   g$oeta <- 0

   OBL <- NULL
   ETA <- NULL
   R2 <- NULL
   SD <- NULL
   debug1 <- ifelse(verbose==2, 1, 0)
   for (filt in filterset) {
      ETAi <- NULL
      R2i <- NULL
      SDi <- NULL
      for (obl in oblset) {
         g$obl <- obl
         model0 <- C_peakFinder(spec, ppmrange, g$flist[[filt]], g, verbose = 0)
         if (model0$nbpeak==0) {
            R2i <- c( R2i, 0 ); SDi <- c( SDi, 1e999 )
            next
         }
         R2j <- NULL
         SDj <- NULL
         for (eta in etaset) {
            g$obl <- obl
            g$peaks <- model0$peaks
            g$peaks$eta <- eta
            model <- C_peakOptimize(spec, ppmrange, g, verbose = 0)
            if (model$nbpeak>0) {
               Ymodel <- model$model + intern_computeBL(spec, model)
               residus <- spec$int-Ymodel
               R2j <- c( R2j, stats::cor(spec$int[iseq],Ymodel[iseq])^2 )
               SDj <- c( SDj, stats::sd(residus[iseq]/spec$Noise) )
            } else {
               R2j <- c( R2j, 0 ); SDj <- c( SDj, 1e999 )
            }
if (debug1) cat(filt,": R2j =",round(R2j[length(R2j)],4),", SDj =",round(SDj[length(SDj)],4)," OBL =",obl," ETA =",eta,"\n")
         }
         idx <- ifelse ( g$criterion==0, which(R2j==max(R2j))[1], which(SDj==min(SDj))[1] )
         ETAi <- c( ETAi, etaset[idx] )
         R2i <- c( R2i, R2j[idx] )
         SDi <- c( SDi, SDj[idx] )
if (debug1) cat(filt,": R2i =",round(R2j[idx],4),", SDi =",round(SDj[idx],4)," OBL =",obl," ETA =",etaset[idx],"\n")
      }
      idx <- ifelse ( g$criterion==0, which(R2i==max(R2i))[1], which(SDi==min(SDi))[1] )
      OBL <- c( OBL, oblset[idx] )
      ETA <- c( ETA, ETAi[idx] )
      R2 <- c( R2, R2i[idx] )
      SD <- c( SD, SDi[idx] )
if (debug1) cat(filt,": R2 =",round(R2i[idx],4),", SD =",round(SDi[idx],4)," OBL =",oblset[idx]," ETA =",ETAi[idx],"\n\n")
   }
   gc()

   R2[ is.na(R2) ] <- 0; SD[ is.na(SD) ] <- 1e999
   if (is.null(R2) || sum(R2)==0)
      stop("No peak found.")

   idx <- ifelse ( g$criterion==0, which(R2==max(R2)), which(SD==min(SD)) )
   fidx <- filterset[idx]

   if (debug1) cat("Best: idx =",idx,", filter =",fidx,", obl =",OBL[idx],"\n")

   model0 <- C_peakFinder(spec, ppmrange, g$flist[[fidx]], g, verbose = debug1)
   g$peaks <- model0$peaks
   g$obl <- OBL[idx];
   g$eta <- ETA[idx];
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
   model$filter <- fidx
   model$crit <- g$crit

   if (verbose) {
      cat('FacN =',FacN,', RatioPN =',g$ratioPN,', RatioSN =',g$ratioPN*FacN,"\n")
      cat('crit =',model$crit,', filter =', model$filter,', obl =',model$params$obl,', eta =',g$eta,"\n")
      cat('Nb Blocks =',model$blocks$cnt,', Nb Peaks =', model$nbpeak,"\n")
      cat('R2 =', model$R2,"\n")
      cat('Residue/Noise : SD =',round(sd(model$residus[iseq])/spec$Noise,4),
                      ', Mean =',round(mean(model$residus[iseq])/spec$Noise,4), "\n")
   }
   class(model) = "LSDmodel"
   model
}

# Local Spectra Deconvolution with predefined peaks (g$peaks not NULL)
LSDeconv_2 <- function(spec, ppmrange, params=NULL, oblset=1:12, verbose=1)
{
   g <- getDeconvParams(params)
   iseq <- getIndexSeq(spec,ppmrange)
   if (is.null(g$facN))
      FacN <- max(20,min(100,round(max(spec$int[iseq])*0.05/spec$Noise)))
   else
      FacN <- g$facN
   spec$B <- spec$Noise/FacN
   g$ratioSN <- FacN*g$ratioPN

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

   if (verbose) {
      cat('FacN =',FacN,', RatioPN =',g$ratioPN,', RatioSN =',g$ratioPN*FacN,"\n")
      cat('crit =',model$crit,', obl =',model$params$obl,"\n")
      cat('Nb Blocks =',model$blocks$cnt,', Nb Peaks =', model$nbpeak,"\n")
      cat('R2 =', model$R2,"\n")
      cat('Residue/Noise : SD =',round(sd(model$residus[iseq])/spec$Noise,4),
                      ', Mean =',round(mean(model$residus[iseq])/spec$Noise,4), "\n")
   }
   class(model) = "LSDmodel"
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
   ppm <- model$ppmrange
   cmodel <- model$peaks[model$peaks$amp/spec$Noise > SNthreshold, -5 ]
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
#' @param title title of the graphic
plotSpec <- function(ppmrange, x, y, ynames=c('Origin', 'Filtered', 'Model'), 
                     ycolors=c('grey', 'blue', 'red', 'green', 'orange','magenta','cyan','darkgreen', 'darkorange'), title='')
{
   iseq <- c(which(x>=ppmrange[1])[1]:length(which(x<=ppmrange[2])))
   if ("numeric" %in% class(y)) {
      data <- data.frame(x=x[iseq], y=y[iseq])
      p <- plotly::plot_ly(data, x = ~x, y = ~y, name = ynames[1], type = 'scatter', mode = 'lines')
   }
   else {
      y1 <- y[,1]
      data <- data.frame(x=x[iseq], y=y1[iseq])
      p <- plotly::plot_ly(data, x = ~x, y = ~y, name = ynames[1], type = 'scatter', mode = 'lines')
      for (k in 2:dim(y)[2]) {
         df <- data.frame(x=x[iseq], y=y[iseq,k])
         p <- p %>% plotly::add_trace(data=df, x = ~x, y = ~y, name=ynames[k], mode = 'lines')
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
   npk <- dim(P2)[1]
   npk_colors <- sample(grDevices::rainbow(npk, s=0.8, v=0.75))
   V <- simplify2array(lapply(1:npk, function(i) { PVoigt(spec$ppm[iseq], P2$amp[i], P2$ppm[i], P2$sigma[i], P2$sigma2[i], P2$eta[i]) }))
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
     data <- data.frame(lab=rownames(model$peaks), x=model$peaks$ppm, y=1.05*model$peaks$amp)
     p1 <- p1 %>% add_annotations(x = data$x, y = data$y, text = data$lab, showarrow = TRUE, arrowcolor='red', 
                           font = list(color = 'black', family = 'sans serif', size = 18))
   }
   p1
}
