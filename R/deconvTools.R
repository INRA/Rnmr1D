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
fsavgol5 <- list( type=3, m=5, nl=5, nr=10 )
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
#'   \item \code{facN} : Peak/Noise Ratio factor - default value = NULL
#'   \item \code{ratioPN} : Peak/Noise Ratio - default value = 1
#'   \item \code{obl} : Optimization of a baseline (BL) for each massif. 0 means no BL, an integer greater than 0 indicates the polynomial order of the BL default value = 0
#'   \item \code{dist_fac} : PeakFinder : min distance between 2 peaks (as multiple of sigma_min which is typically equal to 0.0005 ppm) - default value = 2
#'   \item \code{optim} : Indicates if optimisation is applied - default value = 1
#'   \item \code{oppm} : Indicates if ppm optimisation is applied - default value = 1
#'   \item \code{osigma} : Indicates if sigma optimisation is applied - default value = 1
#'   \item \code{d2spmeth} : PeakFinder : Indicates if minima method to the second derivation is applied
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
  ratioPN = 1,

  # Optimization of a baseline (BL) for each massif
  # 0 means no BL, an integer greater than 0 indicates the polynomial order of the BL
  obl = 0,

  # Peaks searching : min distance between 2 peaks (as multiple of sigma_min which is typically equal to 0.0005 ppm)
  dist_fac = 2,

  # Optimization of peaks : 0 => No, 1 => Yes
  oppm = 1,
  osigma = 1,

  # Peaks searching : Minima method applied to the second derivation
  d2spmeth = 1,

  # CV on Spectrum and its second derivate
  spcv = 0.005,
  d2cv = 0.05,

  # Apply Filter (1) on 1st (d1filt) and 2nd derivates (d2filt) or not (0)
  d1filt = 0,
  d2filt = 1,

  # Optimization of Sigmas : Fixe the limits (min and max) of sigmas
  sigma_min = 0.0005,
  sigma_max = 0.005,

  # Indicates if we want information messages
  verbose = 1,

  # Exclude ppm zones for the criterion evaluation
  exclude_zones = NULL
)

# Merge specific parameters values with the full deconvolution list

#' getDeconvParams
#'
#' \code{getDeconvParams} merges some specific parameters values with the full deconvolution list and return the resulting list. With no parameter as input, it returns the default parameter list.
#' @param params a list defining some specific parameters for deconvolution
#' @return the resulting list of deconvolution parameters.
getDeconvParams <- function(params=NULL)
{
  fullparams <- deconvParams
  if ("list" %in% class(params)) for (p in ls(params)) fullparams[[p]] <- params[[p]]
  fullparams
}

# get the index sequence corresponding to the ppm range
getseq <- function(spec, ppm)
{
   c(which(spec$ppm>=ppm[1])[1]:length(which(spec$ppm<=ppm[2])))
}

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
   out <- C_fSavGol(s, m, nl, nr)
   out
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
   out <- C_FilterbyThreshold(s, wavelet, threshold = type)
   out
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
   out <- C_FilterbyWT(s, type = wavelet, threshold)
   out
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
   out <- C_Lorentz(ppm, amp, x0, sigma)
   out
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
   out <- C_OneLorentz(X, Y, par)
   out
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
#' @param peaks a list of specific parameters for deconvolution, including the matrix defining peaks, one peak by row, with columns defined as : pos, ppm, amp, sigma, pfac, integral
#' @param verbose level of debug information
#' @return a list
peakOptimize <- function(spec, ppmrange, peaks, verbose=1)
{
   model <- C_peakOptimize(spec, ppmrange, peaks, verbose)
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
   out <- C_specModel(spec, ppmrange, peaks)
   out
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
   model <- peakOptimize(spec, ppmrange, params, verbose = verbose)
   class(model) = "optimModel"
   model
}

# Estimated number of peaks
estimation_nbpeaks <- function(spec, ppmrange, params=NULL)
{
   g <- getDeconvParams(params)
   FacN <- ifelse(is.null(g$facN), 5, g$facN)
   spec$B <- spec$Noise/FacN
   Fpars <- list(ratioSN=FacN*g$RationPN, spcv=g$spcv, d2cv=g$d2cv, d2spmeth=g$d2spmeth,
                 d1filt=g$d1filt, d2filt=g$d2filt, dist_fac=g$dist_fac)
   model0 <- peakFinder(spec, ppmrange, Fpars, 'daub8', verbose = 0)
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
#' @return spec object
readSpectrum <- function(ACQDIR, procParams, ppmnoise=c(10.2,10.5), PHC=NULL, scaleIntensity=1)
{
   if (!is.null(PHC)) {
      procParams$OPTPHC0 <<- FALSE
      procParams$OPTPHC1 <<- FALSE
      procParams$phc0 <<- PHC[1]*pi/180
      procParams$phc1 <<- PHC[2]*pi/180
   }
   spec <- Spec1rDoProc(Input=ACQDIR,param=procParams)
   
   Noise <- C_noise_estimation(spec$int, which(spec$ppm>=ppmnoise[1])[1], length(which(spec$ppm<=ppmnoise[2])))
   spec$B <- spec$Noise <- Noise/scaleIntensity
   spec$LAYER_MAX <- round(log(1024*round(length(spec$int)/1024))/log(2)+.05)
   spec$int <- spec$int/scaleIntensity
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
estimateBL <- function(spec, ppmrange, WS=50, NEIGH=35) {
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
   FacN <- ifelse( is.null(g$facN), 5, g$facN )
   spec$B <- spec$Noise/FacN
   debug1 <- ifelse(verbose==2, 1, 0)

   # Peak search
   Fpars <- list(ratioSN=FacN*g$ratioPN, spcv=g$spcv, d2cv=g$d2cv, d2spmeth=g$d2spmeth,
                 d1filt=g$d1filt, d2filt=g$d2filt, dist_fac=g$dist_fac)
   model0 <- peakFinder(spec, ppmrange, Fpars, filter, verbose = debug1)

   # Peak optimization
   Opars <- list(optim=1, os=g$osigma, op=g$oppm, obl=0, ratioSN=FacN*g$ratioPN, oneblk=0,
                 sigma_max=g$sigma_max, tol=g$reltol)
   Opars$peaks <- model0$peaks

   for (k in scset) {
      if (debug1) cat(k,':')
      Opars$scmin <- k
      if (k==min(scset)) {
         model1 <- peakOptimize(spec, ppmrange, Opars, verbose = debug1)
      } else {
         model2 <- peakOptimize(spec, ppmrange, Opars, verbose = debug1)
         model1$peaks <- getBestPeaks(spec, model1, model2, crit=g$crit)
      }
   }
   Ymodel <- specModel(spec, ppmrange, model1$peaks)
   iseq <- getseq(spec,ppmrange)

   model1$model <- Ymodel
   model1$R2 <- stats::cor(spec$int[iseq],Ymodel[iseq])^2
   model1$residus <- spec$int-Ymodel
   model1$SD <- stats::sd(model1$residus[iseq]/spec$Noise)

   if (verbose) {
      cat('FacN =',FacN,', RatioPN =',g$ratioPN,', RatioSN =',g$ratioPN*FacN,"\n")
      cat("Nb Blocks =",model1$blocks$cnt,",Nb Peaks =", model1$nbpeak,"\n")
      cat("R2 =", model1$R2,"\n")
      cat("SD/N =",model1$SD, "\n")
   }
   class(model1) = "GSDmodel"
   model1
}

#=====================================================================
# LSD - Local Spectra Deconvolution
#=====================================================================

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

#' LSDeconv
#'
#' Local Spectra Deconvolution: \code{LSDeconv} belongs to the low-level functions group for deconvolution.
#' @param spec a 'spec' object
#' @param ppmrange a ppm range as a list in order to apply the deconvolution
#' @param filterset a set of filter type for filtering the noise and  smoothing the signal
#' @param oblset a set of baseline order for fitting
#' @param params a list of specific parameters for deconvolution
#' @param verbose level of debug information
#' @return a model object
LSDeconv <- function(spec, ppmrange, params=NULL, filterset=1:6, oblset=1:12, verbose=1)
{
   g <- getDeconvParams(params)
   if (is.null(g$facN))
      FacN <- max(20,min(100,round(max(spec$int[getseq(spec,ppmrange)])*0.05/spec$Noise)))
   else
      FacN <- g$facN
   spec$B <- spec$Noise/FacN

   Opars <- list(ratioSN=FacN*g$ratioPN, spcv=g$spcv, d2cv=g$d2cv, d2spmeth=g$d2spmeth,
                 d1filt=g$d1filt, d2filt=g$d2filt, dist_fac=g$dist_fac, estimate_int=0,
                 sigma_max=g$sigma_max, os=g$osigma, op=g$oppm, obl=1, oneblk=1, tol=g$reltol)

   R2 <- NULL
   SD <- NULL
   OBL <- NULL
   debug1 <- ifelse(verbose==2, 1, 0)
   iseq <- getseq(spec,ppmrange)
   for (filter in filterset) {
      R2i <- NULL
      SDi <- NULL
      for (obl in oblset) {
         Opars$obl <- obl
         modelF <- specDeconv(spec, ppmrange, Opars, filter, verbose = 0)
         lb <- computeBL(spec, modelF)
         Ymodel <- modelF$model + lb
         residus <- spec$int-Ymodel
         R2i <- c( R2i, stats::cor(spec$int[iseq],Ymodel[iseq])^2 )
         SDi <- c( SDi, stats::sd(residus[iseq]/spec$Noise) )
      }
      idx <- ifelse ( g$crit==0, which(R2i==max(R2i)), which(SDi==min(SDi)) )
      OBL <- c( OBL, oblset[idx] )
      R2 <- c( R2, R2i[idx] )
      SD <- c( SD, SDi[idx] )
      if (debug1) cat(filter,": R2 =",R2i[idx],", SD =",SDi[idx]," OBL =",oblset[idx],"\n")
   }
   gc()

   idx <- ifelse ( g$crit==0, which(R2==max(R2)), which(SD==min(SD)) )
   fidx <- filterset[idx]

   if (debug1) cat("Best: idx =",idx,", filter =",fidx,", obl =",Opars$obl,"\n")

   model0 <- peakFinder(spec, ppmrange, Opars, fidx, verbose = 0)
   Opars$obl <- OBL[idx]
   Opars$peaks <- model0$peaks
   modelF <- peakOptimize(spec, ppmrange, Opars, verbose = 0)

   modelF$LB <- computeBL(spec, modelF)
   Ymodel <- modelF$model + modelF$LB
   modelF$residus <- spec$int-Ymodel
   modelF$iseq <- iseq
   modelF$ppmrange <- ppmrange
   modelF$R2 <- stats::cor(spec$int[iseq],Ymodel[iseq])^2
   modelF$CV <- stats::sd(modelF$residus[iseq])/mean(modelF$model[iseq])
   modelF$filter <- fidx
   modelF$crit <- g$crit

   if (verbose) {
      cat('FacN =',FacN,', RatioPN =',g$ratioPN,', RatioSN =',g$ratioPN*FacN,"\n")
      cat("crit =",modelF$crit,", filter =", modelF$filter,", obl =",modelF$params$obl,"\n")
      cat("Nb Peaks =", modelF$nbpeak,"\n")
      cat("R2 =", modelF$R2,"\n")
      cat("CV =",modelF$CV, "\n")
   }
   class(modelF) = "LSDmodel"
   modelF
}


#' cleanPeaks
#'
#' \code{cleanPeaks} cleans the peaks under a specified threshold 
#' @param spec a 'spec' object
#' @param model a model object
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
#' @param spec a 'spec' object (see \code{readSpectrum}, \code{Spec1rDoProc})
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
#' \code{plotModel} plots the model along with the resulting lorentzians from deconvolution
#' @param spec a 'spec' object (see \code{readSpectrum}, \code{Spec1rDoProc})
#' @param model a 'model' object (see \code{specDeconv}, \code{peakOptimize}, \code{GSDeconv}, \code{LSDeconv})
#' @param exclude_zones a list of vector defining the excluded zones for lorentzian plots
#' @param title title of the graphic
plotModel <- function(spec, model, exclude_zones=NULL, title='')
{
   if ( ! sum(c('peakModel', 'optimModel','GSDmodel', 'LSDmodel') %in% class(model) ) )
      stop("the input model must have an appropriate class, namely one of these: 'peakModel', 'optimModel', 'GSDmodel', 'LSDmodel'")
   ppmrange <- c(model$params$wmin, model$params$wmax)
   iseq <- getseq(spec,ppmrange)
   P1 <- model$peaks[model$peaks$ppm>ppmrange[1], ]
   P2 <- P1[P1$ppm<ppmrange[2],]
   if (! is.null(exclude_zones))
     for (k in 1:length(exclude_zones)) P2 <- rbind( P2[P2[,2]<exclude_zones[[k]][1],], P2[P2[,2]>exclude_zones[[k]][2],] )

   npk <- dim(P2)[1]
   npk_colors <- sample(grDevices::rainbow(npk, s=0.8, v=0.75))
   V <- simplify2array(lapply(1:npk, function(i) { Lorentz(spec$ppm[iseq], P2$amp[i], P2$ppm[i], P2$sigma[i]) }))
   fmodel <- apply(V,1,sum)
   datamodel <- data.frame(x=spec$ppm[iseq], ymodel=fmodel)
   p1 <- plotly::plot_ly(datamodel, x = ~x, y = ~ymodel, name = 'Model', type = 'scatter', mode = 'lines' )
   for (i in 1:npk){
     df <- data.frame(x=datamodel$x, y=V[,i])
     p1 <- plotly::add_trace(p1, data=df, x = ~x, y = ~y, name=paste0("p",P2[i,2]), mode = 'lines', fill = 'tozeroy')
   }
   p1 <- p1 %>% plotly::layout(title = title, xaxis = list(autorange = "reversed"), colorway = c('black', npk_colors))
   p1
}
