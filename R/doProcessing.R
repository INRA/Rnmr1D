# ID doProcessing.R
# Copyright (C) 2017-2019 INRA
# Authors: D. Jacob
#

#' detectCores
#'
#' \code{detectCores} is simply a shortcut for parallel::detectCores().
#'
#' @param ...  See \code{parallel::detectCores}
detectCores <- function(...) {
   parallel::detectCores(...)
}

#' setLogFile
#'
#' \code{setLogFile} allows to redirect all log messages to a file
#'
#' @param con  a connection object which inherits from class "connection" 
#' @examples
#'  \donttest{
#'   # Redirect all log messages to a temporary file
#'    outtmp <- tempfile()
#'    con <- file(outtmp, "wt", encoding = "UTF-8")
#'    setLogFile(con)
#'    data_dir <- system.file("extra", package = "Rnmr1D")
#'    RAWDIR <- file.path(data_dir, "CD_BBI_16P02")
#'    CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
#'    SAMPLEFILE <- file.path(data_dir, "Samples.txt")
#'    out <- Rnmr1D::doProcessing(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=6)
#'    close(con)
#'    readLines(outtmp)
#'  }
setLogFile <- function(con=stdout())
{
   if ("connection" %in% class(con) ) LOGFILE <<- con
}

#' setPPMbounds
#'
#' Set the PPM bounds for proton (1H) and carbon (13C) to consider in the processing step and then to store in the specMat object
#'
#' @param proton Minimal and Maximal ppm value for 1H NMR
#' @param carbon Minimal and Maximal ppm value for 13C NMR
setPPMbounds <- function(proton=c(-0.5,11), carbon=c(0,200))
{
   PPM_MIN <<- proton[1]
   PPM_MAX <<- proton[2]
   PPM_MIN_13C <<- carbon[1]
   PPM_MAX_13C <<- carbon[2]
}
   
#' doProcessing 
#'
#' \code{doProcessing} is the main function of this package. Indeed, this function performs 
#' the complete processing of a set of 1D NMR spectra from the FID (raw data) and based on a 
#' processing sequence (macro-command file). An additional file specifies all the spectra to 
#' be considered by associating their sample code as well as the levels of experimental 
#' factors to which they belong. In this way it is possible to select only a subset of spectra 
#' instead of the whole set. 
#'
#' @param path  The full path of either the raw spectra directory on the disk
#' @param cmdfile The full path name of the Macro-commands file for processing (text format)
#' @param samplefile The full path name of the Sample file (tabular format)
#' @param bucketfile The full path name of the file of bucket's zones (tabular format)
#' @param ncpu The number of cores [default: 1]
#' @return
#' \code{doProcessing} returns a list containing the following components:
#' \itemize{
#'   \item \code{samples} : the samples matrix with the correspondence of the raw spectra, 
#' as well as the levels of the experimental factors if specified in the input.
#'   \item \code{factors} : the factors matrix with the corresponding factor names. 
#' At minimum, the list contains the Samplecode label corresponding to the samples without their 
#' group level.
#'   \item \code{rawids} : list of the full directories of the raw spectra (i.e. where the FID files 
#' are accessible)
#'   \item \code{infos} : list of the acquisition and processing parameters for each (raw) spectra.
#'   \item \code{specMat} : objects list  regarding the spectra data.
#'       \itemize{
#'             \item \code{int} : the matrix of the spectra data (\code{nspec} rows X \code{size} 
#' columns)
#'             \item \code{nspec} : the number of spectra
#'             \item \code{size} : the size (i.e number of points) of each spectra
#'             \item \code{ppm_min}, \code{ppm_max} : the minimum and the maximum ppm values of 
#' spectra
#'             \item \code{ppm} : the vector of the ppm values (\code{size} values)
#'             \item \code{dppm} : the ppm increment between each point
#'             \item \code{buckets_zones} : the matrix of the buckets zones including two columns 
#' (min and max) 
#'             \item \code{namesASintMax} : boolean - If TRUE, generate all output matrix with 
#' bucket names based on ppm values of the maximum of the average intensity of all spectra within
#' the ppm range of each bucket. If FALSE (default), then bucket names will be based on the ppm 
#' range center of each bucket.
#'         }
#' }
#' @examples
#'  \donttest{
#'     data_dir <- system.file("extra", package = "Rnmr1D")
#'     cmdfile <- file.path(data_dir, "NP_macro_cmd.txt")
#'     samplefile <- file.path(data_dir, "Samples.txt")
#'     out <- Rnmr1D::doProcessing(data_dir, cmdfile=cmdfile, 
#'                                 samplefile=samplefile, ncpu=detectCores())
#' }
#' @seealso the NMRProcFlow online documentation \url{https://nmrprocflow.org/} and especially 
#' the Macro-command Reference Guide (\url{https://nmrprocflow.org/themes/pdf/Macrocommand.pdf})
#'
doProcessing <- function (path, cmdfile, samplefile=NULL, bucketfile=NULL, ncpu=1 )
{
   if( ! file.exists(path))
       stop(paste0("ERROR: ",path," does NOT exist\n"), call.=FALSE)
   if( ! file.exists(cmdfile))
       stop(paste0("ERROR: ",cmdfile," does NOT exist\n"), call.=FALSE)
   if( ! is.null(samplefile) && ! file.exists(samplefile))
       stop(paste0("ERROR: ",samplefile," does NOT exist\n"), call.=FALSE)
   if( ! is.null(bucketfile) && ! file.exists(bucketfile))
       stop(paste0("ERROR: ",bucketfile," does NOT exist\n"), call.=FALSE)
   if ( checkMacroCmdFile(cmdfile) == 0 )
       stop(paste0("ERROR: ",cmdfile," seems to include errors\n"), call.=FALSE)

   trim <- function (x) gsub("^\\s+|\\s+$", "", x)

   Write.LOG(LOGFILE, "Rnmr1D:  --- READING and CONVERTING ---\n", mode="at")

   # Initialize the list of processing parameters
   procParams <- Spec1rProcpar
   procParams$VENDOR <- "bruker"
   procParams$INPUT_SIGNAL <- "1r"
   procParams$READ_RAW_ONLY <- TRUE

   # Rnmr1D macrocommand file: Get the preprocessing parameter line if exists
   CMDTEXT <- gsub("\t", "", readLines(cmdfile))
   if ( length(grep("#%%", CMDTEXT[1]))==1 ) {
        procpar <- unlist(strsplit(gsub("#%% ", "", CMDTEXT[1]), "; "))
        Write.LOG(LOGFILE, paste0( "Rnmr1D:  ", paste(procpar,collapse=", "), "\n"))
        parnames <- NULL; parvals <- NULL
        for (param in procpar ) { parnames <- c( parnames, unlist(strsplit(param,"="))[1] ); parvals <- c( parvals, unlist(strsplit(param,"="))[2] ); }
        names(parvals) <- parnames;  procpar <- data.frame(t(parvals), stringsAsFactors=FALSE)
        procParams$READ_RAW_ONLY <- FALSE
        if (! is.null(procpar$Vendor)) procParams$VENDOR <- tolower(trim(procpar$Vendor))
        if (! is.null(procpar$Type)) procParams$INPUT_SIGNAL <- trim(procpar$Type)
        if (! is.null(procpar$LB)) procParams$LB <- as.numeric(procpar$LB)
        if (! is.null(procpar$GB)) procParams$GB <- as.numeric(procpar$GB)
        if (! is.null(procpar$BLPHC)) procParams$BLPHC <- ifelse( procpar$BLPHC=="TRUE", TRUE, FALSE)
        if (! is.null(procpar$ZF)) procParams$ZEROFILLING <- as.numeric(procpar$ZF)
        if (! is.null(procpar$PHC0)) procParams$OPTPHC0 <- ifelse( procpar$PHC0=="TRUE", TRUE, FALSE)
        if (! is.null(procpar$PHC1)) procParams$OPTPHC1 <- ifelse( procpar$PHC1=="TRUE", TRUE, FALSE)
        if (! is.null(procpar$ZNEG)) procParams$RABOT <- ifelse( procpar$ZNEG=="TRUE", TRUE, FALSE)
        if (! is.null(procpar$TSP)) procParams$TSP <- ifelse( procpar$TSP=="TRUE", TRUE, FALSE)
        if (! is.null(procpar$O1RATIO)) procParams$O1RATIO <- as.numeric(procpar$O1RATIO)
        if (! is.null(procpar$TSPSNR)) procParams$TSPSNR <- as.numeric(procpar$TSPSNR)
        if (! is.null(procpar$FP)) procParams$FRACPPM <- as.numeric(procpar$FP)
   }

   # Generate the 'samples.csv' & 'factors' files from the list of raw spectra
   Write.LOG(LOGFILE, "Rnmr1D:  Generate the 'samples' & 'factors' files from the list of raw spectra\n")

   metadata <- NULL
   samples <- NULL
   if (!is.null(samplefile) && file.exists(samplefile))
      samples <- utils::read.table(samplefile, sep="\t", header=T,stringsAsFactors=FALSE)

   metadata <- generateMetadata(path, procParams, samples)

   # If ERROR occurs ...
   if (is.null(metadata)) {
       msg <- "Something failed when attempting to generate the metadata files"
       stop(paste0(msg,"\n"), call.=FALSE)
   }
   gc()

   LIST <- metadata$rawids
   Write.LOG(LOGFILE, paste0("Rnmr1D:  -- Nb Spectra = ",dim(LIST)[1]," -- Nb Cores = ",ncpu,"\n"))

   tryCatch({

       cl <- parallel::makeCluster(ncpu)
       doParallel::registerDoParallel(cl)
       Sys.sleep(1)

       x <- 0
       specList <- foreach::foreach(x=1:(dim(LIST)[1]), .combine=cbind) %dopar% {
            ACQDIR <- LIST[x,1]
            NAMEDIR <- ifelse( procParams$VENDOR=='bruker', basename(dirname(ACQDIR)), basename(ACQDIR) )
            PDATA_DIR <- ifelse( procParams$VENDOR=='rs2d', 'Proc', 'pdata' )
            # Init the log filename
            procParams$LOGFILE <- LOGFILE
            procParams$PDATA_DIR <- file.path(PDATA_DIR,LIST[x,3])
            spec <- Spec1rDoProc(Input=ACQDIR,param=procParams)
            if (procParams$INPUT_SIGNAL=='1r') Sys.sleep(0.3)
            Write.LOG(stderr(),".")
            if (dim(LIST)[1]>1) {
                list( x, spec )
            } else {
                spec
            }
       }
       Write.LOG(LOGFILE,"\n")
       gc()

       parallel::stopCluster(cl)

       if (dim(LIST)[1]>1) {
          # Ensure that the specList array is in the same order than both  samples and IDS arrays
          L <- simplify2array(sapply( order(simplify2array(specList[1,])), function(x) { specList[2,x] } ) )
          specList <- L
       }

       Write.LOG(LOGFILE, "Rnmr1D:  Generate the final matrix of spectra...\n")

       M <- NULL
       N <- dim(LIST)[1]
       vpmin<-0; vpmax<-0

       for(i in 1:N) {
           if (N>1) { spec <- specList[,i]; } else { spec <- specList; }
           if (spec$acq$NUC == "13C") { PPM_MIN <- PPM_MIN_13C; PPM_MAX <- PPM_MAX_13C; }
           P <- spec$ppm>PPM_MIN & spec$ppm<=PPM_MAX
           V <- spec$int[P]
           vppm <- spec$ppm[P]
           if (PPM_MIN<spec$pmin) {
               nbzeros <- round((spec$pmin - PPM_MIN)/spec$dppm)
               vpmin <- vpmin + spec$pmin - nbzeros*spec$dppm
               V <- c( rep(0,nbzeros), V )
           } else {
               vpmin <- vpmin + vppm[1]
           }
           if (PPM_MAX>spec$pmax) {
               nbzeros <- round((PPM_MAX - spec$pmax)/spec$dppm)
               vpmax <- vpmax + spec$pmax + nbzeros*spec$dppm
               V <- c( V, rep(0,nbzeros) )
           } else {
               vpmax <- vpmax + vppm[length(vppm)]
           }
           M <- rbind(M, rev(V))
       }

       cur_dir <- getwd()

       specMat <- NULL
       specMat$int <- M
       specMat$ppm_max <- (vpmax/N);
       specMat$ppm_min <- (vpmin/N);
       specMat$nspec <- dim(M)[1]
       specMat$size <- dim(M)[2]
       specMat$dppm <- (specMat$ppm_max - specMat$ppm_min)/(specMat$size - 1)
       specMat$ppm <- rev(seq(from=specMat$ppm_min, to=specMat$ppm_max, by=specMat$dppm))
       specMat$buckets_zones <- NULL
       specMat$namesASintMax <- FALSE  # FALSE -> center, TRUE -> intMax
       specMat$fWriteSpec <- FALSE

       specObj <- metadata
       specObj$specMat <- specMat

       samples <- metadata$samples

       # Raw IDs : expno & procno 
       IDS <- cbind(basename(dirname(as.vector(LIST[,1]))), LIST[, c(2:3)])
       if (N>1) {
             PARS <- t(sapply( c(1:N), function(x) {
                       c( specList[,x]$acq$PULSE, specList[,x]$acq$NUC,   specList[,x]$acq$SOLVENT,    specList[,x]$acq$GRPDLY, 
                          specList[,x]$proc$phc0, specList[,x]$proc$phc1, specList[,x]$acq$SFO1,       specList[,x]$proc$SI, 
                          specList[,x]$acq$SW,    specList[,x]$acq$SWH,   specList[,x]$acq$RELAXDELAY, specList[,x]$acq$O1 )
             }))
             specObj$nuc <- specList[,1]$acq$NUC
       } else {
             PARS <- t( c( spec$acq$PULSE, spec$acq$NUC, spec$acq$SOLVENT, spec$acq$GRPDLY, spec$proc$phc0,      spec$proc$phc1,
                           spec$acq$SFO1,  spec$proc$SI, specList$acq$SW,  spec$acq$SWH,    spec$acq$RELAXDELAY, spec$acq$O1) )
             specObj$nuc <- spec$acq$NUC
       }
       LABELS <- c("PULSE", "NUC", "SOLVENT", "GRPDLY", "PHC0","PHC1","SF","SI","SW", "SWH","RELAXDELAY","O1" )
       
       if (regexpr('BRUKER', toupper(specList[,1]$acq$INSTRUMENT))>0) {
          if (N>1) { PARS <- cbind( IDS[,c(2:3)], PARS ) } else { PARS <- c( IDS[1,c(2:3)], PARS ) }
          LABELS <- c("EXPNO", "PROCNO", LABELS)
       }
       if (regexpr('RS2D', toupper(specList[,1]$acq$INSTRUMENT))>0) {
          if (N>1) { PARS <- cbind( IDS[,3], PARS ) } else { PARS <- c( IDS[1,3], PARS ) }
          LABELS <- c("PROCNO", LABELS)
       }
       if (N>1) {
          PARS <- cbind( samples[,1], samples[,2], PARS )
       } else {
          PARS <- c( samples[1, 1], samples[1, 2], PARS )
       }
       colnames(PARS) <- c("Spectrum", "Samplecode", LABELS )
       specObj$infos <- PARS
       specObj$origin <- paste(procParams$VENDOR, procParams$INPUT_SIGNAL)

     # Rnmr1D processing macrocommand file
       Write.LOG(LOGFILE,"Rnmr1D: ------------------------------------\n")
       Write.LOG(LOGFILE,"Rnmr1D: Process the Macro-commands file\n")
       Write.LOG(LOGFILE,"Rnmr1D: ------------------------------------\n")
       Write.LOG(LOGFILE,"Rnmr1D: \n")

     # Process the Macro-commands file
       specMat <- doProcCmd(specObj, CMDTEXT, ncpu=ncpu, debug=TRUE)
       if (specMat$fWriteSpec) specObj$specMat <- specMat
       gc()

     # Performs the bucketing based on the bucket list file
       if( ! is.null(bucketfile) && file.exists(bucketfile)) {
           Write.LOG(LOGFILE, "Rnmr1D: ------------------------------------\n")
           Write.LOG(LOGFILE, "Rnmr1D: Process the file of buckets\n")
           Write.LOG(LOGFILE, "Rnmr1D: ------------------------------------\n")
           Write.LOG(LOGFILE, "Rnmr1D: \n")

           # The bucket zones file
           buckets_infile <- utils::read.table(bucketfile, header=T, sep="\t",stringsAsFactors=FALSE)
           if ( sum(c('min','max') %in% colnames(buckets_infile)) == 2 ) {
                buckets <- cbind( buckets_infile$max, buckets_infile$min )
                specObj$specMat$buckets_zones <- buckets
                Write.LOG(LOGFILE, paste0("Rnmr1D:     NB Buckets = ",dim(buckets)[1],"\n"))
                Write.LOG(LOGFILE, "Rnmr1D: \n")
           } else {
                Write.LOG(LOGFILE,"ERROR: the file of bucket's areas does not contain the 2 mandatory columns having 'min' and 'max' in its header line\n")
           }
       }

       specObj$specMat$fWriteSpec <- NULL
       specObj$specMat$LOGMSG <- NULL

   }, error=function(e) {
       cat(paste0("ERROR: ",e))
   })

   return(specObj)

}
