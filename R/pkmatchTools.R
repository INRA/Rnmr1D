# ID pkmatchTools.R
# Copyright (C) 2015-2026 INRAE
# Authors: D. Jacob, M. Lefebvre
#

## Advice: Better if implemented by using the R object oriented  programming
## http://www.cyclismo.org/tutorial/R/s4Classes.html

# Parameter List
PM_pars <<- list (
     dataREF  = NULL,   # A reference peaklist for all compound in JSON format
     peaklist = NULL,   # A peak list in JSON format, i.e.  {"peaklists": [{"ppm": value, "i": value}, { }, ..., { } ]}
     delta    = NULL,   # The chemical shift tolerance
     method   = NULL    # Retrieve candidates with at least one peak ("one") or with all peaks ("all")
)

# ---------------------------------------------------------
# Get Parameter
# ---------------------------------------------------------
PM_get <- function(ParName)
{
   return (PM_pars[[ParName]]);
}


# ---------------------------------------------------------
# Transform each char value in numeric value (except ID)
# - return the same list with numerical value
# ---------------------------------------------------------

.doNum <- function(liste) {
    for (i in 1:length(liste$spectraList)) {
        liste$spectraList[[i]]$peaklists[[1]] <- rapply(liste$spectraList[[i]]$peaklists[[1]], as.numeric, how="replace")
    }
    return(liste)
}

# ---------------------------------------------------------
# Return ID + ppm + intensity associated
# ---------------------------------------------------------

.findID <- function (liste, id) {
    toSent <- list()
    for (i in 1:length(liste)) {
        if (liste[[i]][1]$id == id) {
            ppm <- c()
            for (k in seq(1,length(liste[[i]]$peaklists),1)){
                ppm<- c(ppm, (liste[[i]]$peaklists[[k]]$ppm))
            }
            int <- c()
            for (l in seq(1,length(liste[[i]]$peaklists),1)){
                int<- c(int, (liste[[i]]$peaklists[[l]]$i))
            }
            toSent[[id]] <- list("ppm"=sort(ppm), "i"=int[order(ppm)])
        }
    }
    if (length(toSent)==0) {
        stop("No ID found")
    }else {
        return(toSent)
    }
}

# ---------------------------------------------------------
# MATCH CANDIDATES
# - return : a list of candidates for each peak of the peaklist
# ---------------------------------------------------------
PM_findPPM <- function ()
{
    listeREF <- PM_pars$dataREF$spectraList
    peaklist <- PM_pars$peaklist
    delta <- PM_pars$delta
    toSend <- list()
    listIDtoPPM <- list()
    results <- list()
    for (i in 1:length(listeREF)) {
        listIDtoPPM <- c(listIDtoPPM,.findID(listeREF, listeREF[[i]]$id))
    }
    for (e in 1:length(peaklist)) { # for each peak of the peaklist
        id <- c()
        ppm <- c()
        distance <-c()
        pkquery <-c()
        intensity <- c()
        for (j in 1:length(listIDtoPPM)) { # for each compound ID
            for (k in 1:length(listIDtoPPM[[j]]$ppm)) { # for each peak of the ppm list
                if (!is.na(listIDtoPPM[[j]]$ppm[k]) && !is.null(listIDtoPPM[[j]]$ppm[k])) {
                    if ((listIDtoPPM[[j]]$ppm[k] > peaklist[[e]]$ppm-delta) && (listIDtoPPM[[j]]$ppm[k] < peaklist[[e]]$ppm+delta)) {
                        distance  <- c( distance,  round(abs(peaklist[[e]]$ppm-listIDtoPPM[[j]]$ppm[k]),4) )
                        pkquery   <- c( pkquery,   peaklist[[e]]$ppm )
                        id        <- c( id,        names(listIDtoPPM[j]))
                        ppm       <- c( ppm,       listIDtoPPM[[j]]$ppm[k])
                        intensity <- c( intensity, listIDtoPPM[[j]]$i[k])
                    }
                }
            }
        }
        # results[e] <- list(toSend)
        toSend[e] <- list(c(distance=list(distance), pkquery=list(pkquery), id=list(id), ppm=list(ppm), intensity=list(intensity)))
    }
    # toSend <- data.frame(id, ppm, distance, pkquery)
    return(toSend)
}

# ---------------------------------------------------------
# SORT CANDIDATES
# - sort match peak by compound id
# - return list
# ---------------------------------------------------------
PM_SortedByID <- function (candidatesList)
{
    results <- list()
    interm <- list()
    if (PM_pars$method == "one") {
        for (i in 1:length(candidatesList)) {
            for (j in 1:length(candidatesList[[i]]$id)) {
                if (is.null(candidatesList[[i]]$id)) next
                interm[candidatesList[[i]]$id[j]] <- list(c(
                                                    list(c(
                                                        "ppm"=candidatesList[[i]]$ppm[j],
                                                        "dist"=candidatesList[[i]]$distance[j],
                                                        "pkref"=candidatesList[[i]]$pkquery[j],
                                                        "i"=candidatesList[[i]]$intensity[j]
                                                    )), interm[[candidatesList[[i]]$id[j]]]
                                                ))
            }
        }
        results = c(interm,results)
    } else {
        temp <- candidatesList[[1]]$id
        if (length(candidatesList) > 1) {
            for (i in 2:length(candidatesList)) {
                temp <- intersect( candidatesList[[i]]$id, temp )
            }
        }
        #if (! identical(temp, character(0))) {
        if (length(temp)>0) {
            for (i in 1:length(candidatesList)) {
                for (j in 1:length(candidatesList[[i]]$id)) {
                   if (is.null(candidatesList[[i]]$id)) next
                    for (k in 1:length(temp)) {
                        if (temp[k] == candidatesList[[i]]$id[j]) {
                            interm[candidatesList[[i]]$id[j]] <- list(c(
                                                                list(c(
                                                                    "ppm"=candidatesList[[i]]$ppm[j],
                                                                    "dist"=candidatesList[[i]]$distance[j],
                                                                    "pkref"=candidatesList[[i]]$pkquery[j],
                                                                    "i"=candidatesList[[i]]$intensity[j]
                                                                )), interm[[candidatesList[[i]]$id[j]]]
                                                            ))
                        }
                    }
                }
            }
        }else {
            results <- 0
            return(results)
            stop
        }
        results = c(interm,results)
    }
    return(results)
}


# ---------------------------------------------------------
# SORT CANDIDATES DISTANCE
# - sort distance for each compound id
# - return list
# ---------------------------------------------------------

PM_SortCandidates <- function (candidate, id)
{
    peaklist <- PM_pars$peaklist
    delta <- PM_pars$delta
    cand_ppm <- c()
    cand_dist <- c()
    cand_pkref <- c()
    cand_int <- c()
    results <- list()
    for(duo in 1:length(candidate)) {
        cand_ppm <- c(candidate[[duo]][ names(candidate[[duo]]) == "ppm" ][[1]], cand_ppm)
        cand_dist <- c(candidate[[duo]][ names(candidate[[duo]]) == "dist" ][[1]], cand_dist)
        cand_pkref <- c(candidate[[duo]][ names(candidate[[duo]]) == "pkref" ][[1]], cand_pkref)
        cand_int <- c(candidate[[duo]][ names(candidate[[duo]]) == "i" ][[1]], cand_int)
    }
    V <- unique(sort(cand_pkref))
    ppm <- c(); distance <- c(); pkquery <- c(); intensity <- c()
    for (i in 1:length(V)) {
       pkrank <- cand_pkref %in% cand_pkref[ cand_pkref == V[i] ]
       distrank <- order(cand_dist[ pkrank ])
       ppm <- c( ppm, (cand_ppm[ pkrank ])[ distrank ])
       distance <- c( distance, (cand_dist[ pkrank ])[ distrank ])
       intensity <- c( intensity, (cand_int[ pkrank ])[ distrank ])
       pkquery <- c( pkquery, cand_pkref[ pkrank ])
    }

    # Compute the score
    score <- ( length(unique(pkquery)) - (mean(as.numeric(distance))/delta) )/( 1+length(peaklist) )
    score <- round(score, digits=4)

    results[id] <- list(c(ppm=list(sort(ppm)), distance=list(distance[order(ppm)]), pkquery=list(pkquery[order(ppm)]), intensity=list(intensity[order(ppm)]), score=list(score)))
    return(results)
}

# ---------------------------------------------------------
# RETRIEVE ALL REF PPM AND REF INTENSITY LIST FOR CANDIDATES ONLY
# - take reference library and candidates (results of ppm peak matching)
# - return ordered list
# ---------------------------------------------------------
PM_RefByID <- function(refList, candidates, minD=0.04)
{
    sortRef <- list()
    for (i in 1:length(refList$spectraList)) {
        for (j in 1:length(names(candidates))) {
            ref <- refList$spectraList[[i]]
            nameCand <- names(candidates)[j]
            if (ref$id == nameCand) {
                ppm <- c()
                int <- c()
                queryMin <- min(candidates[[j]]$ppm)
                queryMax <- max(candidates[[j]]$ppm)
                for (y in 1:length(ref$peaklists)) {
                    if ( queryMin-minD < as.numeric(ref$peaklists[[y]]$ppm) &&
                        as.numeric(ref$peaklists[[y]]$ppm) < queryMax+minD )
                    {
                        ppm <- c(as.numeric(ref$peaklists[[y]]$ppm), ppm)
                        int <- c(as.numeric(ref$peaklists[[y]]$i), int)
                    }
                }
                sortRef[nameCand] <- list(c( ppm=list(sort(ppm)), i=list(int[order(ppm)]) ))
            }

        }
    }
    return(sortRef)
}

# ---------------------------------------------------------
# SORT CANDIDATES SCORE (Nmatch / 1+Ntotal)
# - return ordered list
# ---------------------------------------------------------
PM_sortScore <- function(candidates)
{
    sortId <- names(candidates)
    scores <- c()
    results <- list()
    for (i in 1:length(candidates)) {
        scores <- c(scores, candidates[[i]]$score)
    }
    #sortID <- sortId[order(scores, decreasing = FALSE)]
    sortID <- sortId[order(scores, decreasing = TRUE)]
    for (j in 1:length(sortID)){
        results[[sortID[j]]] <- candidates[[sortID[j]]]
    }
    return(results)
}

PM_process <- function(dataREF, peaklist, delta, method)
{
    results <- NULL
    PM_pars <<- list (
         dataREF  = .doNum(dataREF),
         peaklist = peaklist,
         delta    = delta,
         method   = method
    )
    repeat {
        candidates <- PM_findPPM()
        sorted_id <- PM_SortedByID(candidates)
        if (class(sorted_id) != "list" || length(sorted_id)==0) break
        sorted_candidates <- list()
        for (i in 1:length(sorted_id))
             sorted_candidates[names(sorted_id[i])] <- PM_SortCandidates(sorted_id[[i]], names(sorted_id[i]))
        results <- PM_sortScore(sorted_candidates)
        break
    }
    return(results)
}

# ---------------------------------------------------------
# Extracts peaklists from the DB compound list (DBNAMES)
# - return a spectraList
# ---------------------------------------------------------
getSpectraList <- function(DBNAMES)
{
    spectraList <- list()
    for (k in 1:nrow(DBNAMES)) {
        db <- as.character(DBNAMES[k,])
        pklist <- unlist(strsplit(as.character(db[5]),","))
        peaklists <- list()
        for(j in 1:length(pklist))
            peaklists[[j]]<- list("ppm"=as.numeric(pklist[j]), "i"="0")
        spectraList[[k]] <- list("id"=db[1], "name"=db[2], "dbid"=db[3], "dblink"=db[4], "peaklists"=peaklists)
    }
     list(spectraList=spectraList)
}

getRefDBinfo <- function(DBNAMES, dbrefid)
{
    dbinfo <- ""
    repeat {
      if (! (dbrefid %in% DBNAMES[,1]) ) break
      dbinfo <- as.character(DBNAMES[which(DBNAMES[,1] == dbrefid),])
      break
  }
  return(dbinfo)
}


#' PeakMatching
#'
#' Searching for putative compounds in a peaklist database from a list of chemical shifts, 
#' obtained for example from a clustering obtained using the \code{getClusters} function.
#'
#' @param peaklist a list of ppm, each corresponding to the maximal intensity of a compound peak 
#' @param lib either an internal library name (dbref6 or hmdb110), or the path of a DB file
#' @param ppmtol ppm tolerance for seaching peaks
#' @param method  either 'all' for all query peaks must machted, or 'one' for only one peak must matched
PeakMatching <- function(peaklist, lib="dbref6",  ppmtol=0.02, method="all")
{
    results <- NULL

    DBfile <- lib
    if (lib %in% c('dbref6', 'hmdb110'))
        DBfile <- paste0(file.path(system.file("extra", package = "Rnmr1D"), toupper(lib)),'_NMR_peaklist.txt')
    if (! file.exists(DBfile))
        stop(paste0(DBfile," does not exist !!"))

    DBNAMES <- read.table(DBfile, sep="\t", header=T)

    dataREF <- getSpectraList(DBNAMES)

    pklist <- list()
    for (i in 1:length(peaklist)) pklist[[i]] <- list("ppm"=peaklist[i])

    results <- PM_process(dataREF, pklist, ppmtol, method)

    for (i in 1:length(results)) {
         dbinfo <- getRefDBinfo(DBNAMES, names(results[i]));
         results[[i]]$name <- dbinfo[2]
         results[[i]]$dbid <- dbinfo[3]
         results[[i]]$dblink <- dbinfo[4]
    }
    results
}
