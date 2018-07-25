#' generateMetadata
#'
#' \code{generateMetadata} Generate the metadata from the list of raw spectra namely the samples, the experimental factors and the list of selected raw spectra. Depending on whether the sample matrix is supplied as input or not, 
#'
#' @param RAWDIR  The full path of either the raw spectra directory on the disk
#' @param procParams the list of processing parameters. First initialize this list with the \code{Spec1r.Procpar.default} list, then modify parameters depending of your spectra set.
#' @param samples the samples matrix with the correspondence of the raw spectra
#'
#' @return
#' \itemize{
#'   \item \code{samples} : the samples matrix with the correspondence of the raw spectra, as well as the levels of the experimental factors if specified in the input.
#'   \item \code{factors} : the factors matrix with the corresponding factor names. At minimum, the list contains the Samplecode label corresponding to the samples without their group level.
#'   \item \code{rawids} : list of the full directories of the raw spectra (i.e. where the FID files are accessible)
#'}
#'
#' @examples
#' data_dir <- system.file("extra", package = "Rnmr1D")
#' samplefile <- file.path(data_dir, "Samples.txt")
#' samples <- read.table(samplefile, sep="\t", header=TRUE,stringsAsFactors=FALSE)
#' metadata <- generateMetadata(data_dir, procParams=Spec1rProcpar, samples)
#'
generateMetadata <- function(RAWDIR, procParams, samples=NULL)
{
   metadata <- list()
   repeat {
      # Bruker, Varian or nmrML with a provided sample file 
      if (! is.null(samples) ) {
          metadata <- set_Metadata(RAWDIR, procParams, samples )
          break
      }
      # Bruker without sample file 
      if ( procParams$VENDOR=="bruker" ) {
          if ( procParams$INPUT_SIGNAL=="fid" ) {
              metadata <- generate_Metadata_fid(RAWDIR, procParams )
              break
          }
          if ( procParams$INPUT_SIGNAL=="1r" )  {
              metadata <- generate_Metadata_1r(RAWDIR, procParams )
              break
          }
          break
      }
      # Varian or nmrML without sample file 
      else {
          metadata <- set_Metadata(RAWDIR, procParams, samples )
          break
      }
      break
   }
   return(metadata)
}

generate_Metadata_fid <- function(RAWDIR, procParams)
{
   metadata <- list()
   ERRORLIST <- c()
   OKRAW <- 1
   lstfac <- matrix(c(1,"Samplecode"), nrow=1)
   RAWPATH <- gsub("//", "/", RAWDIR)
   LIST <- gsub("//", "/", list.files(path = RAWPATH, pattern = "fid$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE))
   if ( class(LIST)=="character" && length(LIST)==0 ) return(NULL)
   L <- simplify2array(strsplit(LIST,'/'))
   if (class(L) != "matrix") {
      L <- simplify2array(lapply(L, length))
      LM <- round(mean(L),0)
      RMLIST <- c()
      if (sum(L<LM)>0) {
          ERRORLIST <- c( ERRORLIST, LIST[ which(L < LM) ] )
          RMLIST <- c( RMLIST, which(L < LM) )
      }
      if (sum(L>LM)>0) {
          ERRORLIST <- c( ERRORLIST, LIST[ which(L > LM) ] )
          RMLIST <- c( RMLIST, which(L > LM) )
      }
      if (length(RMLIST)>0) LIST <- LIST[ -RMLIST ]
   }
   LIST <- as.data.frame(t(simplify2array(strsplit(LIST,'/'))))

   nDir <- dim(simplify2array(strsplit(RAWPATH,'/')))[1]
   LIST <- LIST[, c(-1:-nDir)]
   nc <- dim(LIST)[2]
   nr <- dim(LIST)[1]
   if (nc<3) {
      LIST <- cbind(LIST[,1],LIST)
   }

   SL <- NULL
   if (nc>3) { SL <- LIST[, c(1:(nc-3)) ]; LIST <- LIST[ , c((nc-2):nc) ]; }
   nc=3

   if (length(levels(LIST[,1]))<nr && length(levels(LIST[,1]))>1 &&
       length(levels(LIST[,2]))<nr && length(levels(LIST[,2]))>1) {
      L <- levels(LIST[,1])
      LIST2 <- NULL
      for (i in 1:length(L)) {
         L2 <- LIST[ LIST[,1]==L[i], ]
         LIST2 <- rbind( LIST2, L2[ L2[,2]==levels(L2[,2])[1], ] )
      }
      LIST <- LIST2
   }
   nr <- dim(LIST)[1]
   MS <- as.matrix(LIST)

   if( !is.null(SL)) { LIST2 <- cbind(SL[c(1: dim(LIST)[1])], LIST); } else { LIST2 <- LIST; }
   nc <- dim(LIST2)[2]
   rawdir <- cbind( sapply(1:nr, function(x){ do.call( paste, c( RAWPATH, as.list(LIST2[x,c(1:(nc-1))]), sep="/")) }), MS[, 2], rep(0,nr) )

   if (length(levels(LIST[,1]))==nr) {
      M <-  MS[, c(1,1) ]
   } else {
      M <-  MS[, c(1,2) ]
   }
   if (nr==1 && class(M)=="character") M <- as.matrix(t(M))

   metadata$ERRORLIST <- ERRORLIST
   if (OKRAW==1) {
      metadata$samples <- M
      metadata$rawids <- gsub("//", "/", rawdir)
      metadata$factors <- lstfac
   }
   return(metadata)
}

generate_Metadata_1r <- function(RAWDIR, procParams)
{
   metadata <- list()
   ERRORLIST <- c()
   OKRAW <- 1
   lstfac <- matrix(c(1,"Samplecode"), nrow=1)
   RAWPATH <- gsub("//", "/", RAWDIR)

   LIST <- gsub("//", "/", list.files(path = RAWPATH, pattern = "1r$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE))
   if ( class(LIST)=="character" && length(LIST)==0 ) return(NULL)
   L <- simplify2array(strsplit(LIST,'/'))
   if (class(L) != "matrix") {
      L <- simplify2array(lapply(L, length))
      LM <- round(mean(L),0)
      RMLIST <- c()
      if (sum(L<LM)>0) {
          ERRORLIST <- c( ERRORLIST, LIST[ which(L < LM) ] )
          RMLIST <- c( RMLIST, which(L < LM) )
      }
      if (sum(L>LM)>0) {
          ERRORLIST <- c( ERRORLIST, LIST[ which(L > LM) ] )
          RMLIST <- c( RMLIST, which(L > LM) )
      }
      if (length(RMLIST)>0) LIST <- LIST[ -RMLIST ]
   }
   LIST <- as.data.frame(t(simplify2array(strsplit(LIST,'/'))))

   # Check if we have a Bruker directory structure
   nDir <- dim(simplify2array(strsplit(RAWPATH,'/')))[1]
   LIST <- LIST[, c(-1:-nDir)]
   nc <- dim(LIST)[2]
   nr <- dim(LIST)[1]
   if (nc<5) {
      LIST <- cbind(LIST[,1],LIST)
   }

   SL <- NULL
   if (nc>5) { SL <- LIST[, c(1:(nc-5)) ]; LIST <- LIST[ , c((nc-4):nc) ]; }
   nc=5

   L <- levels(LIST[,4])
   if (length(L)>1) {
      LIST <- LIST[ LIST[,4]==L[1], ]
   }
   nr <- dim(LIST)[1]

   if (length(levels(LIST[,1]))<nr && length(levels(LIST[,1]))>1 &&
       length(levels(LIST[,2]))<nr && length(levels(LIST[,2]))>1) {
      L <- levels(LIST[,1])
      LIST2 <- NULL
      for (i in 1:length(L)) {
         L2 <- LIST[ LIST[,1]==L[i], ]
         LIST2 <- rbind( LIST2, L2[ L2[,2]==levels(L2[,2])[1], ] )
      }
      LIST <- LIST2
   }
   nr <- dim(LIST)[1]
   MS <- as.matrix(LIST)

   if( !is.null(SL)) { LIST2 <- cbind(SL[c(1: dim(LIST)[1])], LIST); } else { LIST2 <- LIST; }
   nc <- dim(LIST2)[2]
   rawdir <- cbind( sapply(1:nr, function(x){ do.call( paste, c( RAWPATH, as.list(LIST2[x,c(1:(nc-3))]), sep="/")) }), MS[, 2], MS[, 4] )

   if (length(levels(LIST[,1]))==nr) {
      M <-  MS[, c(1,1) ]
   } else {
      M <-  MS[, c(1,2) ]
   }
   if (nr==1 && class(M)=="character") M <- as.matrix(t(M))

   metadata$ERRORLIST <- ERRORLIST
   if (OKRAW==1) {
      metadata$samples <- M
      metadata$rawids <- gsub("//", "/", rawdir)
      metadata$factors <- lstfac
   }
   return(metadata)
}

set_Metadata <- function(RAWDIR, procParams, samples)
{
   if (procParams$VENDOR == "bruker") return (.set_Metadata_Bruker(RAWDIR, procParams, samples))
   if (procParams$VENDOR == "varian") return (.set_Metadata_Varian(RAWDIR, procParams, samples))
   if (procParams$VENDOR == "nmrml")  return (.set_Metadata_nmrML(RAWDIR, procParams, samples))
   if (procParams$VENDOR == "jeol")   return (.set_Metadata_Jeol(RAWDIR, procParams, samples))
   return(NULL)
}

.set_Metadata_Bruker <- function(RAWDIR, procParams, samples)
{
   metadata <- list()
   if (!is.null(samples)) {
      samplesize <- dim(samples)
      nraw <- samplesize[1]
      nbcol <- samplesize[2]
      
      lstfac <- matrix(c(1,"Samplecode"), nrow=1)
      rawdir <- NULL
      ERRORLIST <- c()
      OKRAW <- 1
      
      if (procParams$INPUT_SIGNAL == "1r") {
         LIST <- gsub('//', '/', list.files(path = RAWDIR, pattern = "1r$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE))
         if ( class(LIST)=="character" && length(LIST)==0 ) return(0)
         LIST <- as.data.frame(t(simplify2array(strsplit(LIST,'/'))))
         
         # Check if we have a Bruker directory structure
         nDir <- dim(simplify2array(strsplit(RAWDIR,'/')))[1]
         nc <- dim(LIST)[2]
         if ((nc-nDir)>5) {
             RAWDIR <- do.call( paste, c( RAWDIR, as.list(LIST[1,c((nDir+1):(nc-5))]), sep="/"))
         }
         for (i in 1:nraw) {
             if ( file.exists( paste(RAWDIR, samples[i,1],samples[i,3],"pdata",samples[i,4], "1r", sep="/")) ) {
                 rawdir <- rbind( rawdir, c( paste(RAWDIR, samples[i,1], samples[i,3], sep="/"), samples[i,3], samples[1,4] ) )
             } else {
                 ERRORLIST <- c( ERRORLIST, paste(samples[1,1],samples[i,3],"pdata",samples[i,4], "1r", sep="/") )
                 OKRAW <- 0
             }
         }
      } else {
         LIST <- gsub("//", "/", list.files(path = RAWDIR, pattern = "fid$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE))
         LIST <- as.data.frame(t(simplify2array(strsplit(LIST,'/'))))
         
         # Check if we have a Bruker directory structure
         nDir <- dim(simplify2array(strsplit(RAWDIR,'/')))[1]
         nc <- dim(LIST)[2]
         if ((nc-nDir)>3) {
             RAWDIR <- do.call( paste, c( RAWDIR, as.list(LIST[1,c((nDir+1):(nc-3))]), sep="/"))
         }
         for (i in 1:nraw) {
             if ( file.exists( paste(RAWDIR, samples[i,1],samples[i,3],"fid", sep="/")) ) {
                 rawdir <- rbind( rawdir, c( paste(RAWDIR, samples[i,1], samples[i,3], sep="/"), samples[i,3], 0 ) )
             } else {
                 ERRORLIST <- c( ERRORLIST, paste(samples[i,1],samples[i,3],"fid", sep="/") )
                 OKRAW <- 0
             }
         }
      }
      
      if (nbcol>4) {
          M <- samples[,c(-3:-4)]
          lstfac <- rbind( lstfac, cbind( c(2:(nbcol-3)), colnames(samples)[c(-1:-4)] ) )
      } else {
          M <- cbind(samples[,1], samples[,2])
      }
      
      metadata$ERRORLIST <- ERRORLIST
      if (OKRAW==1) {
         metadata$samples <- M
         metadata$rawids <- gsub("//", "/", rawdir)
         metadata$factors <- lstfac
      }
   }
   return(metadata)
}

.set_Metadata_Varian <- function(RAWDIR, procParams, samples)
{
   lstfac <- matrix(c(1,"Samplecode"), nrow=1)
   rawdir <- NULL
   metadata <- list()
   ERRORLIST <- c()
   OKRAW <- 1

   LIST <- gsub('//', '/', list.files(path = RAWDIR, pattern = "fid$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE))
   if ( class(LIST)=="character" && length(LIST)==0 ) return(0)

#   if (!is.null(SampleFile) && file.exists(SampleFile)) {
#       samples <- read.table(SampleFile, sep="\t", header=T,stringsAsFactors=FALSE)
   if (!is.null(samples)) {
       samplesize <- dim(samples)
       nraw <- samplesize[1]
       nbcol <- samplesize[2]

       LIST <- as.data.frame(t(simplify2array(strsplit(LIST,'/'))))
      
       # Check directory structure
       nDir <- dim(simplify2array(strsplit(RAWDIR,'/')))[1]
       nc <- dim(LIST)[2]
       if ((nc-nDir)>2) {
           RAWDIR <- do.call( paste, c( RAWDIR, as.list(LIST[1,c((nDir+1):(nc-2))]), sep="/"))
       }
      for (i in 1:nraw) {
          if ( file.exists( paste(RAWDIR, samples[i,1],"fid", sep="/")) ) {
              rawdir <- rbind( rawdir, c( paste(RAWDIR, samples[i,1], sep="/"), 0, 0 ) )
          } else {
              ERRORLIST <- c( ERRORLIST, paste(samples[i,1],"fid", sep="/") )
              OKRAW <- 0
          }
      }
      if (nbcol==2) {
          M <- cbind(samples[,1], samples[,2])
      }
      if (nbcol>2) {
         M <- samples
         lstfac <- rbind( lstfac, cbind( c(2:(nbcol-1)), colnames(samples)[c(-1:-2)] ) )
      }
   } else {
      rawdir <- cbind( dirname(LIST), rep(0, length(LIST)), rep(0, length(LIST)) )
      M <- cbind( basename(dirname(LIST)), basename(dirname(LIST)) )
   }

   metadata$ERRORLIST <- ERRORLIST
   if (OKRAW==1) {
      metadata$samples <- M
      metadata$rawids <- gsub("//", "/", rawdir)
      metadata$factors <- lstfac
   }
   return(metadata)
}

.set_Metadata_nmrML <- function(RAWDIR, procParams, samples)
{
   return(.set_Metadata_ext(RAWDIR, procParams, samples, ext="nmrML"))
}

.set_Metadata_Jeol <- function(RAWDIR, procParams, samples)
{
   return(.set_Metadata_ext(RAWDIR, procParams, samples, ext="jdf"))
}

.set_Metadata_ext <- function(RAWDIR, procParams, samples, ext="nmrML")
{
   lstfac <- matrix(c(1,"Samplecode"), nrow=1)
   rawdir <- NULL
   metadata <- list()
   ERRORLIST <- c()
   OKRAW <- 1
   pattern <- paste0('.',ext,'$')
   LIST <- gsub('//', '/', list.files(path = RAWDIR, pattern = pattern, all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE))
   if ( class(LIST)=="character" && length(LIST)==0 ) return(0)

   if (!is.null(samples)) {
      samplesize <- dim(samples)
      nraw <- samplesize[1]
      nbcol <- samplesize[2]

      LIST2 <- as.data.frame(t(simplify2array(strsplit(LIST,'/'))))
      nc <- dim(LIST2)[2]
      for (i in 1:nraw) {
          if (sum( samples[i,1] == LIST2[,nc]) == 1) {
              rawdir <- rbind( rawdir,  c( LIST[ which(samples[i,1] == LIST2[,nc]) ], 0, 0)  )
          } else {
              ERRORLIST <- c( ERRORLIST, samples[i,1] )
              OKRAW <- 0
          }
      }
      if (nbcol==2) {
          M <- cbind(samples[,1], samples[,2])
      }
      if (nbcol>2) {
         M <- samples
         lstfac <- rbind( lstfac, cbind( c(2:(nbcol-1)), colnames(samples)[c(-1:-2)] ) )
      }
   } else {
      rawdir <- cbind( LIST, rep(0, length(LIST)), rep(0, length(LIST)) )
      M <- cbind( basename(LIST), gsub(pattern, "", basename(LIST)) )
   }

   metadata$ERRORLIST <- ERRORLIST
   if (OKRAW==1) {
      metadata$samples <- M
      metadata$rawids <- gsub("//", "/", rawdir)
      metadata$factors <- lstfac
   }
   return(metadata)
}