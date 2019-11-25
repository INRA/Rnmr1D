# ID RnmrTools.R
# Copyright (C) 2017-2019 INRA
# Authors: D. Jacob
#

lbALIGN <- 'align'
lbGBASELINE <- 'gbaseline'
lbBASELINE <- 'baseline'
lbQNMRBL <- 'qnmrbline'
lbAIRPLS <- 'airpls'
lbBIN <- 'binning'
lbCALIB <- 'calibration'
lbNORM <- 'normalisation'
lbFILTER <- 'denoising'
lbWARP <- 'warp'
lbCLUPA <- 'clupa'
lbBUCKET <- 'bucket'
lbZERO <- 'zero'
EOL <- 'EOL'

# returns string w/o leading or trailing whitespace
.N <- function(x) { as.numeric(as.vector(x)) }
.C <- function(x) { as.vector(x) }

Write.LOG <- function(logfile=stdout(), ...) cat(sprintf(...), sep='', file=logfile, append=TRUE)

fitdistr <- function(...) {
   MASS::fitdistr(...)
}

#------------------------------
# airPLS
#------------------------------
.WhittakerSmooth <- function(x,w,lambda,differences=1)
{
  x=matrix(x,nrow = 1, ncol=length(x))
  L=length(x)
  E=Matrix::spMatrix(L,L,i=seq(1,L),j=seq(1,L),rep(1,L))
  D=methods::as(diff(E,1,differences),"dgCMatrix")
  W=methods::as(Matrix::spMatrix(L,L,i=seq(1,L),j=seq(1,L),w),"dgCMatrix")
  background=solve((W+lambda*t(D)%*%D),t((w*x)));
  return(as.vector(background))
}
 
.airPLS <- function(x,lambda=100,porder=1, itermax=8)
{
  
  x = as.vector(x)
  m = length(x)
  w = rep(1,m)
  i = 1
  repeat {
     z = .WhittakerSmooth(x,w,lambda,porder)
     d = x-z
     sum_smaller = abs(sum(d[d<0])) 
     if(sum_smaller<0.001*sum(abs(x))||i==itermax) break
     w[d>=0] = 0
     w[d<0] = exp(i*abs(d[d<0])/sum_smaller)
     w[1] = exp(i*max(d[d<0])/sum_smaller)
     w[m] = exp(i*max(d[d<0])/sum_smaller)
     i=i+1
  }
  return(z) 
}


#------------------------------
# Peak detection for spectra
#------------------------------
# Input parameters
#   - X: spectral dataset in matrix format in which each row contains a single sample
#   - nDivRange: size of a single small segment after division of spectra, Default value: 64
#   - baselineThresh: removal of all the peaks with intensity lower than this threshold, Default value: 50000
# Output parameters
#   - peak lists of the spectra
.detectSpecPeaks <- function (X, nDivRange, scales=seq(1, 16, 2), baselineThresh, SNR.Th=-1) 
{
  nFea <- ncol(X)
  nSamp <- nrow(X)
  noiseEsp <- 0.005
  if (SNR.Th < 0) SNR.Th <- max(scales) * 0.05
  i<-0
  pList <- foreach::foreach(i=1:nSamp) %dopar% {
     myPeakRes <- NULL
     mySpec <- X[i, ]
     for (k in 1:length(nDivRange)) {
        divR <- nDivRange[k]
        for (j in 1:(trunc(nFea/divR) - 3)) {
           startR <- (j - 1) * divR + 1
           if (startR >= nFea)  startR <- nFea
           endR <- (j + 3) * divR
           if (endR > nFea) endR <- nFea
           xRange <- mySpec[startR:endR]
           xMean <- mean(xRange)
           xMedian <- stats::median(xRange)
           if ((xMean == xMedian) || abs(xMean - xMedian)/((xMean + xMedian) * 2) < noiseEsp) next
           peakInfo <- MassSpecWavelet::peakDetectionCWT(mySpec[startR:endR], scales = scales, SNR.Th = SNR.Th)
           majorPeakInfo <- peakInfo$majorPeakInfo
           if (length(majorPeakInfo$peakIndex) > 0) myPeakRes <- c(myPeakRes, majorPeakInfo$peakIndex + startR - 1)
        }
     }
     plst <- list(myPeakRes)
     plst[[1]] <- sort(unique(plst[[1]]))
     plst[[1]] <- plst[[1]][which(mySpec[ plst[[1]] ] > baselineThresh)]
     plst[[1]] <- sort(plst[[1]])
     plst
  }
  pList <-  simplify2array(pList)
  return(pList)
}


#------------------------------
# CluPA alignment for multiple spectra
#------------------------------
# Input parameters
#   - X: spectral dataset in the matrix format in which each row contains a single sample
#   - peakList: peak lists of the spectra
#   - refInd: index of the reference spectrum
#   - maxShift:  maximum number of the points for a shift step
# Output parameters
#   - aligned spectra: same format as input
.dohCluster <- function (X, peakList, refInd = 1, maxShift = 50, acceptLostPeak = TRUE) 
{
  refSpec = X[refInd, ]
  tarInd<-0
  Y <- foreach::foreach(tarInd=1:nrow(X), .combine=rbind) %dopar% {
     if (tarInd == refInd) { refSpec; } else {
     tarSpec <- X[tarInd, ]
     myPeakList <- c(peakList[[refInd]], peakList[[tarInd]])
     myPeakLabel <- double(length(myPeakList))
     myPeakLabel[1:length(peakList[[refInd]])] <- 1
     res <- hClustAlign(refSpec, tarSpec, myPeakList, myPeakLabel, 1, length(tarSpec), maxShift = maxShift, acceptLostPeak = TRUE)
     res$tarSpec
  }}
  return(Y)
}

#------------------------------
# Spectra alignment - see https://cran.r-project.org/web/packages/speaq/vignettes/speaq.pdf
#------------------------------
# Input parameters
#   - data: n x p datamatrix
#   - nDivRange: size of a single small segment after division of the whole spectrum, Default value: 64
#   - reference: number of the spectrum reference; if NULL, automatic detection, Default value: NULL
#   - baselineThresh: removal of all the peaks with intensity lower than this threshold, Default value: 50000
# Output parameters
#   - Y: n x p datamatrix
.CluPA <- function(data, reference=reference, nDivRange, scales = seq(1, 16, 2), baselineThresh,  SNR.Th = -1, maxShift=50, DEBUG=FALSE)
{
  LOGMSG <- ""

  ## Peak picking
  if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste("Rnmr1D:     --- Peak detection : nDivRange =",nDivRange,"\n"));
  startTime <- proc.time()
  peakList <- .detectSpecPeaks(X=data, nDivRange=nDivRange, scales=scales, baselineThresh=baselineThresh, SNR.Th = SNR.Th)
  endTime <- proc.time()
  if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste("Rnmr1D:     --- Peak detection time: ",(endTime[3]-startTime[3])," sec\n"));

  ## Reference spectrum determination
  if (reference == 0) {
     #resFindRef<- speaq::findRef(peakList)
     #refInd <- resFindRef$refInd
     V <- bestref(data, optim.crit="WCC")
     refInd  <- V$best.ref

  } else  {
     refInd=reference
  }
  if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste("Rnmr1D:     --- The reference spectrum is: ",refInd,"\n"));

  ## Spectra alignment to the reference
  if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste("Rnmr1D:     --- Spectra alignment to the reference: maxShift =",maxShift,"\n"));
  startTime <- proc.time()
  Y <- .dohCluster(data, peakList=peakList, refInd=refInd, maxShift=maxShift, acceptLostPeak=TRUE)
  endTime <- proc.time()
  if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste("Rnmr1D:     --- Spectra alignment time: ",(endTime[3]-startTime[3])," sec\n"));

  ## Output  
  return(list(M=Y, LOGMSG=LOGMSG))
}

#------------------------------
# Calibration ot the PPM Scale
#------------------------------
RCalib1D <- function(specMat, PPM_NOISE_AREA, zoneref, ppmref)
{
   i1<-length(which(specMat$ppm>max(zoneref)))
   i2<-which(specMat$ppm<=min(zoneref))[1]

   # PPM calibration of each spectrum
   for ( i in 1:specMat$nspec ) {
       i0 <- i1 + which(specMat$int[i, i1:i2]==max(specMat$int[i, i1:i2])) - 1
       ppm0 <- specMat$ppm_max - (i0-1)*specMat$dppm
       dppmref <- ppm0 - ppmref
       decal <- 0
       sig <- C_estime_sd(specMat$int[i, 1:specMat$size],128)
       if (abs(dppmref) > specMat$dppm) {
           decal <- trunc(dppmref/specMat$dppm)
           dppmref <- dppmref - decal*specMat$dppm
       }
       if (abs(dppmref) > 0.5*specMat$dppm) {
           decal <- decal + trunc(2*dppmref/specMat$dppm)
           dppmref <- dppmref - trunc(2*dppmref/specMat$dppm)*specMat$dppm
       }
       if (decal==0) next

       if (decal<0) {
          specMat$int[i, 1:(specMat$size-abs(decal))] <- specMat$int[i,(1+abs(decal)):specMat$size]
          specMat$int[i, (specMat$size-abs(decal)+1):specMat$size] <- stats::rnorm(length((specMat$size-abs(decal)+1):specMat$size), mean=specMat$int[i,specMat$size-abs(decal)-1], sd=sig)
       }
       if (decal>0) {
          specMat$int[i,(1+abs(decal)):specMat$size] <- specMat$int[i, 1:(specMat$size-abs(decal))]
          specMat$int[i, 1:abs(decal)] <- stats::rnorm(length(1:abs(decal)), mean=specMat$int[i,abs(decal)+1], sd=sig)
       }
   }
   return(specMat)
}

#------------------------------
# Normalisation of the Intensities
#------------------------------
RNorm1D <- function(specMat, normmeth, zones)
{
   N <- dim(zones)[1]

   if (normmeth=='CSN') {
      # 1/ Integration of each zone ...
      i <- 0
      SUM <- foreach::foreach(i=1:N, .combine='+') %dopar% {
          i1<-length(which(specMat$ppm>max(zones[i,])))
          i2<-which(specMat$ppm<=min(zones[i,]))[1]
          simplify2array(lapply( 1:specMat$nspec, function(x) { 
                0.5*(specMat$int[x, i1] + specMat$int[x, i2]) + sum(specMat$int[x,(i1+1):(i2-1)])
          }))
      }
      COEFF <- SUM/mean(SUM)
   }
   if (normmeth=='PQN') {
      # 1/ Get spectra values for each zone ...
      i <- 0
      SUBMAT <- foreach::foreach(i=1:N, .combine=cbind) %dopar% {
          i1<-length(which(specMat$ppm>max(zones[i,])))
          i2<-which(specMat$ppm<=min(zones[i,]))[1]
          # .. for each spectrum
          t(simplify2array(lapply( 1:specMat$nspec, function(x) { specMat$int[x,i1:i2] })))
      }
      # Calculate the most probabe quotient
      V <- apply(SUBMAT, 2, stats::median)
      MQ <- t(t(SUBMAT)/V)
      COEFF <- apply(MQ,1, stats::median)
   }

   # 2/ Apply to each spectrum, its corresponding coefficient
   V <- lapply( 1:specMat$nspec, function(x) { specMat$int[x,] <<- specMat$int[x,]/COEFF[x] } )
   return(specMat)
}

#------------------------------
# Global Baseline Correction
#------------------------------
RGbaseline1D <- function(specMat,PPM_NOISE_AREA, zone, WS, NEIGH)
{

   NFAC <- 1.5
   NEGFAC <- 10
   NBPASS <- 2
   
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   TD <- specMat$size

   # Baseline Estimation for each spectrum
   i <- 0
   BLList <- foreach::foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       sig <- fitdistr(specMat$int[i, length(which(specMat$ppm>PPM_NOISE_AREA[2])):(which(specMat$ppm<=PPM_NOISE_AREA[1])[1])], "normal")$estimate[2]
       mv <- simplify2array(lapply( c(3:61), function(x) { mean(abs(specMat$int[i, ((x-1)*TD/64):(x*TD/64)])); }))
       mmoy<-min(mv); mv <- NULL
       if (min(specMat$int[i, ])< -NEGFAC*mmoy) {
           Vthreshold <- max(specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ])
           specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ] <- Vthreshold
       }
       V <- specMat$int[i, ]
       if (NBPASS==1) {
           BL <- C_Estime_LB (V, i1, i2, WS, NEIGH, NFAC*sig)
       } else {
           BL <- 0*rep(1:length(V))
           for( n in 1:NBPASS) {
              # Estimation of Baseline
              BLn <- C_Estime_LB (V, i1, i2, WS, NEIGH, NFAC*sig)
              V <- V - BLn
              BL <- BL + BLn
           }
       }
       BL
   }
   gc()

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       mv <- simplify2array(lapply( c(3:61), function(x) { mean(abs(specMat$int[i, ((x-1)*TD/64):(x*TD/64)])); }))
       mmoy<-min(mv); mv <- NULL
       if (min(specMat$int[i, ])< -NEGFAC*mmoy) {
           Vthreshold <- max(specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ])
           specMat$int[ i, specMat$int[i, ]/mmoy < -NEGFAC ] <- Vthreshold
       }
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       specMat$int[i,c(i1:i2)] <- specMat$int[i,c(i1:i2)] - BL[c(i1:i2)]
   }

   return(specMat)
}

#------------------------------
# Local Baseline Correction (deprecated)
#------------------------------
RBaseline1D <- function(specMat,PPM_NOISE_AREA, zone, WINDOWSIZE)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )

   SMOOTHSIZE <- round(WINDOWSIZE/2)
   ALPHA <- 0.2

   # Ajust the window size parameter
   n <- i2-i1+1
   ws <- WINDOWSIZE
   dws <- 0
   signdws <- ifelse ( n/ws>round(n/ws), 1, -1 )
   dmin <- 1
   if (n>WINDOWSIZE) {
      repeat {
          d <- abs(n/(ws+signdws*dws)-round(n/(ws+signdws*dws)))
          if ( d>dmin ) { dws <- dws - signdws; break }
          dmin <- d; dws <- dws + signdws;
      }
      WINDOWSIZE <- ws+signdws*dws
      n2 <- round(n/(ws+signdws*dws))*(ws+signdws*dws)
      if (n2<n) i2 <- i2 - (n-n2)
   } else {
      WINDOWSIZE <- n
   }

   # Baseline Estimation for each spectrum
   i<-0
   BLList <- foreach::foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       specSig <- fitdistr(specMat$int[i,length(which(specMat$ppm>PPM_NOISE_AREA[2])):(which(specMat$ppm<=PPM_NOISE_AREA[1])[1])], "normal")$estimate[2]
       x <- specMat$int[i,c(i1:i2)]
       xmat <- matrix(x, nrow=WINDOWSIZE)
       ymin <- apply(xmat, 2, min) + 1.2*specSig
       y1 <- rep(ymin, each = WINDOWSIZE)
       y2 <- Smooth(y1,SMOOTHSIZE)
       BL <- simplify2array(lapply(c(1:length(x)), function(k) { min(y1[k], y2[k]); }))
       BL <- lowpass1(BL, ALPHA)
       BL
   }

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       specMat$int[i,c(i1:i2)] <- specMat$int[i,c(i1:i2)] - BL
   }

   return(specMat)
}

#------------------------------
# q-NMR Baseline Correction
#------------------------------
Rqnmrbc1D <- function(specMat, PPM_NOISE_AREA, zone)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   n <- i2-i1+1
   NLOOP <- 5
   dN <- round(0.00075/specMat$dppm)
   CSIG <- 5

   # Baseline Estimation for each spectrum
   i<-0
   BLList <- foreach::foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       specSig <- fitdistr(specMat$int[i,length(which(specMat$ppm>PPM_NOISE_AREA[2])):(which(specMat$ppm<=PPM_NOISE_AREA[1])[1])], "normal")$estimate[2]
       x <- specMat$int[i,c(i1:i2)]
       bc <- 0*rep(1:n)
       for (l in 1:NLOOP) {
           bci <- C_GlobSeg(x, dN, CSIG*specSig)
           x <- x - bci
           bc <- bc + bci
       }
       bc
  }

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       V <- specMat$int[i,c(i1:i2)] - BL
       specMat$int[i,c(i1:i2)] <- V
   }

   return(specMat)
}

#------------------------------
# airPLS : Local Baseline Correction
#------------------------------
RairPLSbc1D <- function(specMat, zone, clambda)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   n <- i2-i1+1
   cmax <- 6

   lambda <- ifelse (clambda==cmax, 5, 10^(cmax-clambda) )
   # Baseline Estimation for each spectrum
   BLList <- foreach::foreach(i=1:specMat$nspec, .combine=cbind) %dopar% {
       x <- specMat$int[i,c(i1:i2)]
       bc <- .airPLS(x, lambda)
       gc()
       bc
   }

   # Baseline Correction for each spectrum
   for ( i in 1:specMat$nspec ) {
       if( is.null(dim(BLList)) ) { BL <- BLList; } else { BL <- BLList[,i]; }
       V <- specMat$int[i,c(i1:i2)] - BL
       specMat$int[i,c(i1:i2)] <- V
   }

   return(specMat)
}

#------------------------------
# Denoising the selected PPM ranges
#------------------------------
RFilter1D <- function(specMat,zone, FILTORD, FILTLEN)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )

   sgfilt <- sgolay(p=FILTORD, n=FILTLEN)

   # Denoising each spectrum
   for ( i in 1:specMat$nspec ) {
        x <- specMat$int[i,c(i1:i2)]
        SpecMat_sg <- filter(sgfilt,x)
        specMat$int[i,c(i1:i2)] <- SpecMat_sg
   }

   return(specMat)

}

#------------------------------
# Zeroing the selected PPM ranges
#------------------------------
RZero1D <- function(specMat, zones, DEBUG=FALSE)
{
   # Zeroing each PPM range
   N <- dim(zones)[1]
   LOGMSG <- ""

   for ( i in 1:N ) {
       i1<-length(which(specMat$ppm>max(zones[i,])))
       i2<-which(specMat$ppm<=min(zones[i,]))[1]
       specMat$int[,c(i1:i2)] <- matrix(0,specMat$nspec,(i2-i1+1))
       if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste("Rnmr1D:     Zone",i,"= (",min(zones[i,]),",",max(zones[i,]),")\n"))
   }
   specMat$LOGMSG <- LOGMSG
   return(specMat)
}

#------------------------------
# LS : Alignment of the selected PPM ranges
#------------------------------
RAlign1D <- function(specMat, zone, RELDECAL=0.35, idxSref=0, Selected=NULL)
{
   # Alignment of each PPM range
   NBPASS <- 3
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   
   decal <- round((i2-i1)*RELDECAL)
   for( n in 1:NBPASS) {
       ret <- align_segment(specMat$int, segment_shifts( specMat$int, idxSref, decal, i1-1, i2-1, Selected-1), i1-1, i2-1, Selected-1)
   }

   return(specMat)
}

#------------------------------
# CluPA : Alignment of the selected PPM ranges
#------------------------------
RCluPA1D <- function(specMat, zonenoise, zone, resolution=0.02, SNR=3, idxSref=0, Selected=NULL, DEBUG=FALSE)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )
   
   # Noise estimation
   PPM_NOISE_AREA <- c(min(zonenoise), max(zonenoise))
   idx_Noise <- c( length(which(specMat$ppm>PPM_NOISE_AREA[2])),(which(specMat$ppm<=PPM_NOISE_AREA[1])[1]) )
   Vref <- spec_ref(specMat$int)
   ynoise <- C_noise_estimation(Vref,idx_Noise[1],idx_Noise[2])
   
   # Parameters
   baselineThresh <- SNR*mean( C_noise_estimate(specMat$int, idx_Noise[1],idx_Noise[2], 1) )
   nDivRange <- max( round(resolution/specMat$dppm,0), 64 )
   maxshift <- min( round(0.01/specMat$dppm), round(nDivRange/4) )

   # Subpart of spectra
   if( is.null(Selected)) M<-specMat$int[, c(i1:i2) ] else  M<-specMat$int[Selected, c(i1:i2) ];

   out <- .CluPA(M, reference=idxSref, nDivRange, scales = seq(1, 8, 2), 
                          baselineThresh,  SNR.Th = 0.1, maxShift=maxshift, DEBUG=DEBUG)

   if( is.null(Selected)) specMat$int[ ,c(i1:i2)] <- out$M else specMat$int[ Selected,c(i1:i2)] <- out$M
   specMat$LOGMSG <- out$LOGMSG

   return(specMat)
}

#------------------------------
# PTW : Alignment of the selected PPM ranges
#------------------------------
RWarp1D <- function(specMat, zone, idxSref=0, warpcrit=c("WCC","RMS"), Selected=NULL)
{
   i1 <- ifelse( max(zone)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone))) )
   i2 <- ifelse( min(zone)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone))[1] )

   if( is.null(Selected)) M<-specMat$int[, c(i1:i2) ] else  M<-specMat$int[Selected, c(i1:i2) ];
   if (is.null(Selected)) nspec <- specMat$nspec else nspec <- length(Selected);

   if (idxSref==0 || ( !is.null(Selected) && !(idxSref %in% Selected) )) {
      V      <- bestref(M, optim.crit=warpcrit)
      refid  <- V$best.ref
   } else {
     refid <- idxSref
   }

   ref    <- M[refid, ]
   ssampl <- M[ c(1:nspec)[-refid], ]
   out    <- ptw(ref, ssampl, warp.type = "individual", mode = "forward", optim.crit=warpcrit)
   #out    <- ptw(ref, ssampl, warp.type = "global", mode = "forward", init.coef = c(0, 1, 0), optim.crit=warpcrit)
   M      <- out$warped.sample
   M[is.na(M)] <- 0
   if( is.null(Selected)) specMat$int[ c(1:nspec)[-refid],c(i1:i2)] <- M else specMat$int[ Selected[-refid],c(i1:i2)] <- M

   return(specMat)
}

#------------------------------
# Bucket : Apply the bucketing based on the 'Algo' algorithm with the resolution 'resol'.
# Then elinate buckets with a SNR under the threshold given by 'snr'
# Append to / or Write upon the bucket file depending the 'appendBuc' value
#------------------------------
RBucket1D <- function(specMat, Algo, resol, snr, zones, zonenoise, appendBuc, DEBUG=FALSE)
{
   # Limit size of buckets
   MAXBUCKETS<-2000
   NOISE_FAC <- 3

   # Noise estimation
   PPM_NOISE_AREA <- c(min(zonenoise), max(zonenoise))
   idx_Noise <- c( length(which(specMat$ppm>PPM_NOISE_AREA[2])),(which(specMat$ppm<=PPM_NOISE_AREA[1])[1]) )
   Vref <- spec_ref(specMat$int)
   ynoise <- C_noise_estimation(Vref,idx_Noise[1],idx_Noise[2])
   Vnoise <- abs( C_noise_estimate(specMat$int, idx_Noise[1],idx_Noise[2], 1) )

   bdata <- list()
   bdata$ynoise <- ynoise
   bdata$vnoise <- NULL
   bdata$inoise_start <- idx_Noise[1]
   bdata$inoise_end <- idx_Noise[2]
   bdata$R <- resol
   bdata$dppm <- specMat$dppm
   bdata$noise_fac <- NOISE_FAC
   bdata$bin_fac <- 0.5
   bdata$peaknoise_rate <- 15
   bdata$BUCMIN <- 0.003
   bdata$VREF <- 1
   LOGMSG <- ""

   # For each PPM range
   buckets_zones <- NULL
   N <- dim(zones)[1]
   i<-0
   buckets_zones <- foreach::foreach(i=1:N, .combine=rbind) %dopar% {
       i2<-which(specMat$ppm<=min(zones[i,]))[1]
       i1<-length(which(specMat$ppm>max(zones[i,])))
       if (Algo=='aibin') {
          Mbuc <- matrix(, nrow = MAXBUCKETS, ncol = 2)
          Mbuc[] <- 0
          buckets_m <- C_aibin_buckets(specMat$int, Mbuc, Vref, bdata, i1, i2)
       }
       if (Algo=='unif') {
          seq_buc <- seq(i1, i2, round(resol/specMat$dppm))
          n_bucs <- length(seq_buc) - 1
          buckets_m <- cbind ( seq_buc[1:n_bucs], seq_buc[2:(n_bucs+1)])
       }
       # Keep only the buckets for which the SNR average is greater than 'snr'
       MaxVals <- C_maxval_buckets (specMat$int, buckets_m)
       buckets_m <- buckets_m[ which( apply(t(MaxVals/(2*Vnoise)),1,stats::quantile)[4,]>snr), ]
       LOGMSG <- paste("Rnmr1D:     Zone",i,"= (",min(zones[i,]),",",max(zones[i,]),"), Nb Buckets =",dim(buckets_m)[1],"\n")
       cbind( specMat$ppm[buckets_m[,1]], specMat$ppm[buckets_m[,2]], LOGMSG )
   }
   if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste(unique(buckets_zones[,3]), collapse=""))

   buckets_zones <- cbind( .N(buckets_zones[,1]), .N(buckets_zones[,2]) )
   if( DEBUG ) LOGMSG <- paste0(LOGMSG, paste("Rnmr1D:     Total Buckets =",dim(buckets_zones)[1],"\n"))

   if( buckets_zones[1,1]>buckets_zones[1,2] )  {  colnames(buckets_zones) <- c('max','min') }
                                          else  {  colnames(buckets_zones) <- c('min','max') }
   if (appendBuc==1) {
      specMat$buckets_zones <- rbind(specMat$buckets_zones, buckets_zones )
   } else {
      specMat$buckets_zones <- buckets_zones
   }
   specMat$LOGMSG <- LOGMSG

   return(specMat)
}

#' checkMacroCmdFile
#'
#' \code{checkMacroCmdFile} Check if the macro-commands included in the input file (commandfile) 
#' are compliant with the allowed commands.
#'
#' @param commandfile  The macro-commands file - the allowed commands are : 'align', 'warp', 
#' 'clupa', 'gbaseline', 'baseline', 'qnmrbline', 'airpls', 'binning', 'calibration', 
#' 'normalisation', 'denoising', 'bucket', 'zero'.
#'
#' @return return 1 if the macro-commands included in the input file are compliant, 0 if not.
#'
#' @examples
#' data_dir <- system.file("extra", package = "Rnmr1D")
#' CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
#' ret <- checkMacroCmdFile(CMDFILE)
#'
#' @seealso the NMRProcFlow online documentation \url{https://nmrprocflow.org/} and especially 
#' the Macro-command Reference Guide (\url{https://nmrprocflow.org/themes/pdf/Macrocommand.pdf})
checkMacroCmdFile <- function(commandfile) {
   ret <- 1
   allowKW <- c( 'align', 'warp', 'clupa', 'gbaseline', 'baseline', 'qnmrbline', 'airpls', 'binning', 'calibration', 'normalisation', 'denoising', 'bucket', 'zero', 'EOL' )

   tryCatch({
      # Read the macrocommand file
      CMDTEXT <- gsub("\t", "", readLines(commandfile))
      CMDTEXT <- CMDTEXT[ grep( "^[^ ]", CMDTEXT ) ]
      CMD <- CMDTEXT[ grep( "^[^#]", CMDTEXT ) ]
      CMD <- gsub("^ ", "", gsub(" $", "", gsub(" +", ";", CMD)))
      L <- unique(sort(gsub(";.*$","", CMD)))
      L <- L[ grep( "^[^0-9]", L)]
      ret <- ifelse( sum(L %in% allowKW)==length(L), 1, 0 )
   }, error=function(e) {
       ret <- 0
   })
   return(ret)
}

#' RWrapperCMD1D
#'
#' \code{RWrapperCMD1D} belongs to the low-level functions group - it serves as a wrapper to 
#' call internale functions for processing.
#'
#' @param cmdName the name of internal function 
#' @param specMat a 'specMat' object
#' @param ... specific parameters of the requested function 
#' @return 
#'  \code{specMat} : a 'specMat' object
RWrapperCMD1D <- function(cmdName, specMat, ...)
{
   repeat {
       if (cmdName == lbCALIB) {
          specMat <- RCalib1D(specMat, ...)
          break
       }
       if (cmdName == lbNORM) {
          specMat <- RNorm1D(specMat, ...)
          break
       }
       if (cmdName == lbGBASELINE) {
          specMat <- RGbaseline1D(specMat, ...)
          break
       }
       if (cmdName == lbBASELINE) {
          specMat <- RBaseline1D(specMat, ...)
          break
       }
       if (cmdName == lbQNMRBL) {
          specMat <- Rqnmrbc1D(specMat, ...)
          break
       }
       if (cmdName == lbAIRPLS) {
          specMat <- RairPLSbc1D(specMat, ...)
          break
       }
       if (cmdName == lbFILTER) {
          specMat <- RFilter1D(specMat, ...)
          break
       }
       if (cmdName == lbZERO) {
          specMat <- RZero1D(specMat, ...)
          break
       }
       if (cmdName == lbALIGN) {
          specMat <- RAlign1D(specMat, ...)
          break
       }
       if (cmdName == lbCLUPA) {
          specMat <- RCluPA1D(specMat, ...)
          break
       }
       if (cmdName == lbWARP) {
          specMat <- RWarp1D(specMat, ...)
          break
       }
       if (cmdName == lbBUCKET) {
          specMat <- RBucket1D(specMat, ...)
          break
       }
       break
    }
   return(specMat)
}

#' doProcCmd
#'
#' \code{doProcCmd} it process the Macro-commands string array specified at input.
#'
#' @param specObj a complex list return by \code{doProcessing} function. See the manual 
#' page of the \code{\link{doProcessing}} function for more details on its structure.
#' @param cmdstr the Macro-commands string array; See the Macro-command Reference Guide 
#' (\url{https://nmrprocflow.org/themes/pdf/Macrocommand.pdf}) to have more details about 
#' macro-commands.
#' @param ncpu The number of cores [default: 1]
#' @param debug a boolean to specify if we want the function to be more verbose.
#' @return 
#'  \code{specMat} : a 'specMat' object - See the manual page of the \code{\link{doProcessing}} 
#' function for more details on its structure
#' @examples
#'  \donttest{
#'     data_dir <- system.file("extra", package = "Rnmr1D")
#'     CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
#'     SAMPLEFILE <- file.path(data_dir, "Samples.txt")
#'     out <- Rnmr1D::doProcessing(data_dir, cmdfile=CMDFILE, 
#'                                 samplefile=SAMPLEFILE, ncpu=detectCores())
#' # Apply an intelligent bucketing (AIBIN)
#'     specMat.new <- Rnmr1D::doProcCmd(out, 
#'              c("bucket aibin 10.2 10.5 0.3 3 0", 
#'                  "9.5 4.9", 
#'                  "4.8 0.5", 
#'                  "EOL"
#'               ),ncpu=2, debug=TRUE)
#'     out$specMat <- specMat.new
#' }
doProcCmd <- function(specObj, cmdstr, ncpu=1, debug=FALSE)
{
 # specMat
   specMat <- specObj$specMat

 # specParams
   specParamsDF <- as.data.frame(specObj$infos, stringsAsFactors=FALSE)
   SI <- as.numeric(specParamsDF$SI[1])

 # Samples
   samples <- specObj$samples

   CMDTEXT <- cmdstr[ grep( "^[^ ]", cmdstr ) ]
   CMD <- CMDTEXT[ grep( "^[^#]", CMDTEXT ) ]
   CMD <- gsub("^ ", "", gsub(" $", "", gsub(" +", ";", CMD)))

   specMat$fWriteSpec <- FALSE

   cl <- parallel::makeCluster(ncpu)
   doParallel::registerDoParallel(cl)
   Sys.sleep(1)

   while ( length(CMD)>0 && CMD[1] != EOL ) {
   
      cmdLine <- CMD[1]
      cmdPars <- unlist(strsplit(cmdLine[1],";"))
      cmdName <- cmdPars[1]

      repeat {
          if (cmdName == lbCALIB) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)>=3) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMREF <- params[3]
                 PPM_NOISE <- ifelse( length(params)==5, c( min(params[4:5]), max(params[4:5]) ), c( 10.2, 10.5 ) )
                 Write.LOG(LOGFILE, paste0("Rnmr1D:  Calibration: PPM REF =",PPMREF,", Zone Ref = (",PPMRANGE[1],",",PPMRANGE[2],")\n"));
                 specMat <- RWrapperCMD1D(cmdName,specMat, PPM_NOISE, PPMRANGE, PPMREF)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbNORM) {
              params <- cmdPars[-1]
              if (length(params)==2) {
                 params <- as.numeric(params)
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Normalisation: Zone Ref = (",PPMRANGE[1],",",PPMRANGE[2],")\n"));
                 specMat <- RWrapperCMD1D(cmdName,specMat, normmeth='CSN', zones=matrix(PPMRANGE,nrow=1, ncol=2))
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              if (length(params)==1) {
                 NORM_METH <- params[1]
                 CMD <- CMD[-1]
                 zones <- NULL
                 while(CMD[1] != EOL) {
                    zones <- rbind(zones, as.numeric(unlist(strsplit(CMD[1],";"))))
                    CMD <- CMD[-1]
                 }
                 Write.LOG(LOGFILE,"Rnmr1D:  Normalisation of the Intensities based on the selected PPM ranges...\n")
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Method =",NORM_METH,"\n"))
                 specMat <- RWrapperCMD1D(cmdName,specMat, normmeth=NORM_METH, zones=zones)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbGBASELINE || cmdName == lbBASELINE) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)==6) {
                 PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMRANGE <- c( min(params[3:4]), max(params[3:4]) )
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Baseline Correction: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )\n"));
                 if (cmdName == lbBASELINE) {
                     if (params[5]>0) { # to be compliant with version<1.1.4
                        BCMETH <- params[5]
                        WSFAC <- (7-params[6])/4
                        WINDOWSIZE <- round(WSFAC*SI/(BCMETH*64))
                     } else {
                        WINDOWSIZE <- round(( 1/2^(params[6]-2) )*(SI/64))
                     }
                     Write.LOG(LOGFILE,paste0("Rnmr1D:     Type=Local - Window Size = ",WINDOWSIZE,"\n"));
                     specMat <- RWrapperCMD1D(cmdName,specMat,PPM_NOISE, PPMRANGE, WINDOWSIZE)
                 } else {
                     WS <- params[5]
                     NEIGH <- params[6]
                     Write.LOG(LOGFILE,paste0("Rnmr1D:     Type=Global - Smoothing Parameter=",WS," - Window Size=",NEIGH,"\n"));
                     specMat <- RWrapperCMD1D(cmdName,specMat,PPM_NOISE, PPMRANGE, WS, NEIGH)
                 }
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbQNMRBL) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)==4) {
                 PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMRANGE <- c( min(params[3:4]), max(params[3:4]) )
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Baseline Correction: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )\n"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Type=q-NMR\n"))
                 specMat <- RWrapperCMD1D(cmdName,specMat,PPM_NOISE, PPMRANGE)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbAIRPLS) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)==3) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 LAMBDA <- params[3]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Baseline Correction: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )\n"))
                 Write.LOG(LOGFILE,paste("Rnmr1D:     Type=airPLS, lambda=",LAMBDA,"\n"))
                 specMat <- RWrapperCMD1D(cmdName,specMat, PPMRANGE, LAMBDA)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }

          if (cmdName == lbFILTER) {
              params <- as.numeric(cmdPars[-1])
              if (length(params)==4) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 FORDER <- params[3]
                 FLENGTH <- params[4]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Denoising: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )\n"));
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Filter Order=",FORDER," - Filter Length=",FLENGTH,"\n"));
                 specMat <- RWrapperCMD1D(cmdName,specMat,PPMRANGE, FORDER, FLENGTH)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbALIGN) {
              params <- as.numeric(cmdPars[-1])
              Selected <- NULL
              if (length(params)==6 && (params[5]<2 || params[6])) {
                 level <- unique(samples[ order(samples[, params[5]+1]), params[5]+1 ])[params[6]]
                 Selected <- .N(rownames(samples[ samples[, params[5]+1]==level, ]))
              }
              if (length(params)>=4) {
                 PPMRANGE <- c( min(params[1:2]), max(params[1:2]) )
                 RELDECAL= params[3]
                 idxSref=params[4]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Alignment: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )\n"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     Rel. Shift Max.=",RELDECAL," - Reference=",idxSref,"\n"))
                 specMat <- RWrapperCMD1D(cmdName,specMat, PPMRANGE, RELDECAL, idxSref, Selected=Selected)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbWARP) {
              params <- cmdPars[-1]
              Selected <- NULL
              if (length(params)==6 && (.N(params[5])<2 || .N(params[6]))) {
                 level <- unique(samples[ order(samples[, .N(params)[5]+1]), .N(params)[5]+1 ])[.N(params[6])]
                 Selected <- .N(rownames(samples[ samples[, .N(params[5])+1]==level, ]))
              }
              if (length(params)>=4) {
                 PPMRANGE <- c( min(.N(params[1:2])), max(.N(params[1:2])) )
                 idxSref=.N(params[3])
                 warpcrit=.C(params[4])
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Alignment: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )\n"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Parametric Time Warping Method - Reference=",idxSref," - Optim. Crit=",warpcrit,"\n"))
                 specMat <- RWrapperCMD1D(cmdName,specMat, PPMRANGE, idxSref, warpcrit, Selected=Selected)
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbCLUPA) {
              params <- as.numeric(cmdPars[-1])
              Selected <- NULL
              if (length(params)==9 && (params[8]<2 || params[9])) {
                 level <- unique(samples[ order(samples[, params[8]+1]), params[8]+1 ])[params[9]]
                 Selected <- .N(rownames(samples[ samples[, params[8]+1]==level, ]))
              }
              if (length(params)>=7) {
                 PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
                 PPMRANGE <- c( min(params[3:4]), max(params[3:4]) )
                 RESOL= params[5]
                 SNR= params[6]
                 idxSref=params[7]
                 Write.LOG(LOGFILE,paste0("Rnmr1D:  Alignment: PPM Range = ( ",min(PPMRANGE)," , ",max(PPMRANGE)," )\n"))
                 Write.LOG(LOGFILE,paste0("Rnmr1D:     CluPA - Resolution =",RESOL," - SNR threshold=",SNR, " - Reference=",idxSref,"\n"))
                 specMat <- RWrapperCMD1D(cmdName,specMat, PPM_NOISE, PPMRANGE, RESOL, SNR, idxSref, Selected=Selected, DEBUG=debug)
                 if (debug) Write.LOG(LOGFILE, specMat$LOGMSG )
                 specMat$fWriteSpec <- TRUE
                 CMD <- CMD[-1]
              }
              break
          }
          if (cmdName == lbZERO) {
              CMD <- CMD[-1]
              zones2 <- NULL
              while(CMD[1] != EOL) {
                  zones2 <- rbind(zones2, as.numeric(unlist(strsplit(CMD[1],";"))))
                  CMD <- CMD[-1]
              }
              Write.LOG(LOGFILE,"Rnmr1D:  Zeroing the selected PPM ranges ...\n")
              specMat <- RWrapperCMD1D(cmdName,specMat, zones2, DEBUG=debug)
              if (debug) Write.LOG(LOGFILE, specMat$LOGMSG )
              specMat$fWriteSpec <- TRUE
              CMD <- CMD[-1]
              break
          }
          if (cmdName == lbBUCKET) {
              if ( length(cmdPars) < 7 || ! cmdPars[2] %in% c('aibin','unif') ) {
                 CMD <- CMD[-1]
                 break;
              }
              params <- as.numeric(cmdPars[-c(1:2)])
              PPM_NOISE <- c( min(params[1:2]), max(params[1:2]) )
              CMD <- CMD[-1]
              zones <- NULL
              while(CMD[1] != EOL) {
                  zones <- rbind(zones, as.numeric(unlist(strsplit(CMD[1],";"))))
                  CMD <- CMD[-1]
              }
              Write.LOG(LOGFILE,"Rnmr1D:  Bucketing the selected PPM ranges ...\n")
              Write.LOG(LOGFILE,paste0("Rnmr1D:     ",toupper(cmdPars[2])," - Resolution =",params[3]," - SNR threshold=",params[4], " - Append=",params[5],"\n"))
              specMat <- RWrapperCMD1D(cmdName,specMat, cmdPars[2], params[3], params[4], zones, PPM_NOISE, params[5], DEBUG=debug)
              specMat$buckets_zones <- specMat$buckets_zones[order(specMat$buckets_zones[,1]), ]
              if (debug) Write.LOG(LOGFILE, specMat$LOGMSG )
              specMat$fWriteSpec <- TRUE
              CMD <- CMD[-1]
              break
          }
          CMD <- CMD[-1]
          break
      }
      gc()
   }

   parallel::stopCluster(cl)

   return(specMat)
}

#' plotSpecMat Overlaid/Stacked Plot
#'
#' \code{plotSpecMat} Plot spectra set, overlaid or stacked; if stacked, plot with or 
#' without a perspective effect.
#'
#' @param specMat a 'specMat' object - Spectra matrix in specMat$int (rows = samples, 
#' columns = buckets)
#' @param ppm_lim ppm range of the plot
#' @param K Graphical height of the stack (0 .. 1),(default=0.67)
#' @param pY Intensity limit factor (default 1)
#' @param dppm_max Max ppm shift to have a perspective effect
#' @param asym Correction of vertical parallax effect  (-1 .. 1)
#'                 -1 : parallelogram
#'                  0 : trapeze with maximum asymmetric
#'                  1 : symmetric trapeze
#' @param beta Correction of horizontal parallax effect   (0 .. 0.2) (defaut 0)
#' @param cols Vector of colors (same size that the number of spectra, i.e dim(specmat)[1])
#'
#' @examples
#'  \donttest{
#'   data_dir <- system.file("extra", package = "Rnmr1D")
#'   cmdfile <- file.path(data_dir, "NP_macro_cmd.txt")
#'   samplefile <- file.path(data_dir, "Samples.txt")
#'   out <- Rnmr1D::doProcessing(data_dir, cmdfile=cmdfile, 
#'                                 samplefile=samplefile, ncpu=detectCores())
#'   # Overlaid plot
#'   plotSpecMat(out$specMat, ppm_lim=c(0.5,9), K=0, pY=0.1)
#'   # Stacked plot with perspective effect
#'   plotSpecMat(out$specMat, ppm_lim=c(-0.1,9),K=0.33)
#'   # Stacked plot with perspective effect with maximum asymmetric
#'   plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33, asym=0)
#'   cols <- c(rep("red",3), rep("blue",3))
#'   # Stacked plot with colors accordings to group levels
#'   plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.67, dppm_max=0, cols=cols)
#' }
plotSpecMat <- function(specMat, ppm_lim=c(min(specMat$ppm),max(specMat$ppm)), K=0.67, pY=1, dppm_max=0.2*(max(ppm_lim) - min(ppm_lim)), asym=1, beta=0, cols=NULL)
{
   specmat <- specMat$int
   ppm <- specMat$ppm
   i2<-which(ppm<=min(ppm_lim))[1]
   i1<-length(which(ppm>max(ppm_lim)))
   ppm_sub <- ppm[i1:i2]
   specmat_sub <- specmat[ , i1:i2 ]
   if (K==0) dppm_max <- 0

   Ymax <- pY*max(specmat_sub)         # intensity limit
   diffppm <- max(ppm_sub) - min(ppm_sub)
   Cy <- 1-beta*(max(ppm_sub)-asym*dppm_max - ppm_sub)/(max(ppm_sub)-asym*dppm_max - min(ppm_sub))
   if (is.null(cols)) cols <- grDevices::rainbow(dim(specmat_sub)[1], s=0.8, v=0.75)

   graphics::plot( cbind(ppm_sub, specmat_sub[1,]), xlim=rev(ppm_lim), ylim=c(0,Ymax), type="h", col="white", xlab="ppm", ylab="Intensity (u.a)")

   graphics::segments( max(ppm_sub) - asym*dppm_max, K*Ymax, min(ppm_sub) + dppm_max, (1-beta)*K*Ymax , col="lightgrey" )
   graphics::segments( max(ppm_sub), 0, min(ppm_sub), 0, col="lightgrey" )
   graphics::segments( max(ppm_sub), 0, max(ppm_sub) - asym*dppm_max, K*Ymax , col="lightgrey" )
   graphics::segments( min(ppm_sub) + dppm_max, (1-beta)*K*Ymax, min(ppm_sub), 0, col="lightgrey" )

   for( i in 1:dim(specmat_sub)[1] ) {
      dppm <- dppm_max*(i-1)/(dim(specmat_sub)[1]-1)
      ppmi <- ppm_sub*(1-(1+asym)*dppm/diffppm) + dppm*(1+(1+asym)*min(ppm_sub)/diffppm)
      y_offset <- K*Ymax*Cy*(i-1)/(dim(specmat_sub)[1]-1)
      graphics::lines(cbind(ppmi, specmat_sub[i,] + y_offset), col=cols[i])
   }
}

#' getBucketsTable
#'
#' Generates the buckets table
#'
#' @param specObj a complex list return by \code{doProcessing} function. See the manual page of the \code{\link{doProcessing}} function for more details on its structure.
#' @return
#'   the buckets table 
getBucketsTable <- function(specObj)
{
   outtable <- NULL
   specMat <- specObj$specMat
   buckets <- specMat$buckets_zones
   if ( ! is.null(buckets) ) {
      colnames(buckets) <- c("max","min")
      buckets_m <- t(simplify2array(lapply( c( 1:(dim(buckets)[1]) ), 
                     function(x) { c( length(which(specMat$ppm>buckets[x,1])), length(which(specMat$ppm>buckets[x,2])) ) }
                    )))
      buckets <- as.data.frame(buckets, stringsAsFactors=FALSE)
      buckets$center <- 0.5*(buckets[,1]+buckets[,2])
      buckets$width <-  0.5*abs(buckets[,1]-buckets[,2])
      buckets$intMax <- specMat$ppm[ C_ppmIntMax_buckets(specMat$int, buckets_m) ]
      if ( is.null(specMat$namesASintMax) || ! specMat$namesASintMax ) {
          buccenter <- buckets$center
      } else {
          buccenter <- buckets$intMax
      }
      buckets$name <- gsub("^(\\d+)","B\\1", gsub("\\.", "_", gsub(" ", "", sprintf("%7.4f",buccenter))) )
      outtable <- buckets[, c("name", "center", "min", "max", "width", "intMax") ]
   }
   return(outtable)
}

#' getBucketsDataset
#'
#' Generates the matrix including the integrations of the areas defined by the buckets (columns) 
#' on each spectrum (rows)
#'
#' Before bucket data export in order to make all spectra comparable with each other, the variations of the overall concentrations of samples have to be taken into account. We propose two normalization methods. In NMR metabolomics, the total intensity normalization (called the Constant Sum Normalization) is often used so that all spectra correspond to the same overall concentration. It simply consists in normalizing the total intensity of each individual spectrum to a same value. An other method called Probabilistic Quotient Normalization (Dieterle et al. 2006) assumes that biologically interesting concentration changes influence only parts of the NMR spectrum, while dilution effects will affect all metabolites signals. Probabilistic Quotient Normalization (PQN) starts by the calculation of a reference spectrum based on the median spectrum. Next, for each variable of interest the quotient of a given test spectrum and reference spectrum is calculated and the median of all quotients is estimated. Finally, all variables of the test spectrum are divided by the median quotient.
#' An internal reference can be used to normalize the data. For example, an electronic reference (ERETIC, see Akoka et al. 1999, or ERETIC2 generated with TopSpin software) can be used for this purpose. The integral value of each bucket will be divided by the integral value of the ppm range given as reference.
#'
#' @param specObj a complex list return by \code{doProcessing} function. See the manual page of the \code{\link{doProcessing}} function for more details on its structure.
#' @param norm_meth Normalization method. The possible values are : 'none', 'CSN' or 'PDN'. See below.
#' @param zoneref Specify the ppm zone of the internal reference (i.e. ERETIC) if applicable. default is NA.
#'
#' @return the data matrix
#'
#' @references{
#'    Akoka S1, Barantin L, Trierweiler M. (1999) Concentration Measurement by Proton NMR 
#' Using the ERETIC Method, Anal. Chem 71(13):2554-7. doi: 10.1021/ac981422i.
#'
#'    Dieterle F., Ross A., Schlotterbeck G. and Senn H. (2006). Probabilistic Quotient 
#' Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. 
#' Application in 1H NMR Metabonomics. Analytical Chemistry, 78:4281-4290.doi: 10.1021/ac051632c
#' }
#'
#' @examples
#'  \donttest{
#'   data_dir <- system.file("extra", package = "Rnmr1D")
#'   cmdfile <- file.path(data_dir, "NP_macro_cmd.txt")
#'   samplefile <- file.path(data_dir, "Samples.txt")
#'   out <- Rnmr1D::doProcessing(data_dir, cmdfile=cmdfile, 
#'                                 samplefile=samplefile, ncpu=detectCores())
#'   outMat <- getBucketsDataset(out, norm_meth='CSN')
#'  }
getBucketsDataset <- function(specObj, norm_meth='none', zoneref=NA)
{
   outdata <- NULL
   specMat <- specObj$specMat
   buckets <- specMat$buckets_zones

   if ( ! is.null(buckets) ) {
      # get index of buckets' ranges
      colnames(buckets) <- c("max","min")
      buckets_m <- t(simplify2array(lapply( c( 1:(dim(buckets)[1]) ), 
                     function(x) { c( length(which(specMat$ppm>buckets[x,1])), length(which(specMat$ppm>buckets[x,2])) ) }
                    )))
      # Integration
      buckets_IntVal <- C_all_buckets_integrate (specMat$int, buckets_m, 0)
      if (norm_meth == 'CSN') {
          buckets_IntVal <- C_buckets_CSN_normalize( buckets_IntVal )
      }
      if (norm_meth == 'PQN') {
          buckets_IntVal_CSN <- C_buckets_CSN_normalize( buckets_IntVal )
          bucVref_IntVal <- C_MedianSpec(buckets_IntVal_CSN)
          bucRatio <- buckets_IntVal_CSN / bucVref_IntVal
          Coeff <- apply(bucRatio,1, stats::median)
          buckets_IntVal <- buckets_IntVal_CSN / Coeff
      }
      # if supplied, integrate of all spectra within the PPM range of the reference signal
      Vref <- 0*c(1:specMat$nspec)
      if (! is.na(zoneref)) {
          istart <- length(which(specMat$ppm>max(zoneref)))
          iend <- length(which(specMat$ppm>min(zoneref)))
          Vref <- C_spectra_integrate (specMat$int, istart, iend)
          buckets_IntVal <- buckets_IntVal/Vref
      }
      # read samples
      samples <- specObj$samples
      # write the data table
      if ( is.null(specMat$namesASintMax) || ! specMat$namesASintMax ) {
          buccenter <- 0.5*(buckets[,1]+buckets[,2])
      } else {
          buccenter <- specMat$ppm[ C_ppmIntMax_buckets(specMat$int, buckets_m) ]
      }
      bucnames <- gsub("^(\\d+)","B\\1", gsub("\\.", "_", gsub(" ", "", sprintf("%7.4f",buccenter))) )
      outdata <- buckets_IntVal
      colnames(outdata) <- bucnames
      rownames(outdata) <- samples[,2]
   }
   return(outdata)
}

#' getSnrDataset
#'
#' Generates the Signal-Noise-Ratio dataset
#'
#' whatever the bucketing approach used, the Signal-to-Noise ratio is a good quality indicator. 
#' Thus, it is possible to check buckets based on their Signal-to-Noise ratio.
#'
#' @param specObj a complex list return by \code{doProcessing} function. 
#' See the manual page of the \code{\link{doProcessing}} function for more details on its structure.
#' @param zone_noise Specify a ppm range of noisy zone default is c(10.2,10.5)
#' @param ratio boolean; TRUE for output Signal-Noise Ratio, or FALSE to output maximum value of 
#' each bucket and in addition, the estimate noise as a separate column
#' @return
#'   the Signal-Noise-Ratio matrix

getSnrDataset <- function(specObj, zone_noise=c(10.2,10.5), ratio=TRUE)
{
   outdata <- NULL
   specMat <- specObj$specMat
   buckets <- specMat$buckets_zones

   if ( ! is.null(buckets) ) {
      # get index of buckets' ranges
      colnames(buckets) <- c("max","min")
      buckets_m <- t(simplify2array(lapply( c( 1:(dim(buckets)[1]) ), 
                     function(x) { c( length(which(specMat$ppm>buckets[x,1])), length(which(specMat$ppm>buckets[x,2])) ) }
                    )))
      # Compute Vnoise vector & Maxvals maxtrix
      i1 <- ifelse( max(zone_noise)>=specMat$ppm_max, 1, length(which(specMat$ppm>max(zone_noise))) )
      i2 <- ifelse( min(zone_noise)<=specMat$ppm_min, specMat$size - 1, which(specMat$ppm<=min(zone_noise))[1] )
      flg <- 1
      Vnoise <- abs( C_noise_estimate(specMat$int, i1, i2, flg) )
      MaxVals <- C_maxval_buckets (specMat$int, buckets_m)

      # read samples
      samples <- specObj$samples
      # write the data table
      if ( is.null(specMat$namesASintMax) || ! specMat$namesASintMax ) {
          buccenter <- 0.5*(buckets[,1]+buckets[,2])
      } else {
          buccenter <- specMat$ppm[ C_ppmIntMax_buckets(specMat$int, buckets_m) ]
      }
      bucnames <- gsub("^(\\d+)","B\\1", gsub("\\.", "_", gsub(" ", "", sprintf("%7.4f",buccenter))) )
      if (ratio) {
         outdata <- cbind( MaxVals/(2*Vnoise))
         colnames(outdata) <- bucnames
      } else {
         outdata <- cbind( Vnoise, MaxVals )
         colnames(outdata) <- c( 'Noise', bucnames )
      }
      rownames(outdata) <- samples[,2]
   }
   return(outdata)
}


#' getSnrDataset
#'
#' Generates the spectral data matrix. The first column indicates the value of ppm, then the following columns correspond to spectral data, one column per spectrum.
#'
#' @param specObj a complex list return by \code{doProcessing} function. 
#' @return
#'   the spectral data matrix
getSpectraData <- function(specObj)
{
   # read samples
   specMat <- specObj$specMat
   samples <- specObj$samples
   outdata <- cbind( specMat$ppm, t(specMat$int) )
   colnames(outdata) <- c( "ppm", samples[,1] )
   return(outdata)
}
