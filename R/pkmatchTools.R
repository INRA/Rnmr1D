# ID pkmatchTools.R
# Copyright (C) 2015-2026 INRAE
# Authors: D. Jacob, M. Lefebvre
#

## Advice: Better if implemented by using the R object oriented  programming
## http://www.cyclismo.org/tutorial/R/s4Classes.html


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
PM_findPPM <- function (PM_pars)
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
PM_SortedByID <- function (PM_pars, candidatesList)
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

PM_SortCandidates <- function (PM_pars, candidate, id)
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
    PM_pars <- list (
         dataREF  = .doNum(dataREF),
         peaklist = peaklist,
         delta    = delta,
         method   = method
    )
    repeat {
        candidates <- PM_findPPM(PM_pars)
        sorted_id <- PM_SortedByID(PM_pars, candidates)
        if (! inherits(sorted_id,"list") || length(sorted_id)==0) break
        sorted_candidates <- list()
        for (i in 1:length(sorted_id))
             sorted_candidates[names(sorted_id[i])] <- PM_SortCandidates(PM_pars, sorted_id[[i]], names(sorted_id[i]))
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

#' PeakMatching
#'
#' Searching for putative compounds in a peaklist database from a list of chemical shifts,
#' based for example on a clustering obtained using the \code{getClusters} function.
#'
#' @param peaklist a list of ppm, each corresponding to the maximal intensity of a compound peak 
#' @param lib either an internal library name (dbref6 or hmdb110), or the path of a DB file
#' @param ppmtol ppm tolerance for seaching peaks
#' @param method  either 'all' for all query peaks must machted, or 'one' for only one peak must matched
peakMatching <- function(peaklist, lib="dbref6",  ppmtol=0.02, method="all")
{
    results <- NULL

    DBfile <- lib
    if (lib %in% c('dbref6', 'hmdb110'))
        DBfile <- paste0(file.path(system.file("extra", package = "Rnmr1D"), toupper(lib)),'_NMR_peaklist.txt')
    if (! file.exists(DBfile))
        stop(paste0(DBfile," does not exist !!"))

    DBNAMES <- utils::read.table(DBfile, sep="\t", header=T)

    dataREF <- getSpectraList(DBNAMES)

    pklist <- list()
    for (i in 1:length(peaklist)) pklist[[i]] <- list("ppm"=peaklist[i])

    results <- PM_process(dataREF, pklist, ppmtol, method)

    for (i in 1:length(results)) {
         dbinfo <- as.character(DBNAMES[which(DBNAMES[,1] == names(results[i])), ])
         results[[i]]$name <- dbinfo[2]
         results[[i]]$dbid <- dbinfo[3]
         results[[i]]$dblink <- dbinfo[4]
         results[[i]]$intensity  <- NULL
    }
    results
}


#' matchClusters
#'
#' Match all clusters obtained using the \code{getClusters} function and using the the \code{peakMatching} function.
#'
#' @param clusters object obtained using the \code{getClusters} function
#' @param lib either an internal library name (dbref6 or hmdb110), or the path of a DB file
#' @param score_min selects only matchs with a score above the value
#' @param ppmtol ppm tolerance for seaching peaks
#' @param method  either 'all' for all query peaks must machted, or 'one' for only one peak must matched
matchClusters <- function(clusters, lib = "dbref6", score_min = 0.5, ppmtol = 0.02, method = "all")
{
    if (! inherits(clusters, "clusters"))
        stop('clusters variable is not a clusters class')
    M <- NULL
    for (k in 1:length(clusters)) {
        C <- paste0("C", k)
        if (!is.null(clusters[[C]])) {
            PK <- NULL
            tryCatch({
                PK <- peakMatching(clusters[[C]], lib,  ppmtol, method)
            }, error=function(e){ PK <- NULL })
            if (!is.null(PK) && length(PK[[1]]) > 0) {
                pk <- PK[[1]]
                if (pk$score > score_min)
                   M <- rbind(M, c(C, pk$name, pk$score, length(unique(pk$pkquery)), length(unique(pk$ppm))))
            }
        }
    }
    colnames(M) <- c("Cluster", "Name", "Score", "Query", "Found")
    out <- as.data.frame(M)

    # Among the groups associated with the same compound, remove those whose score is below the selection threshold (10% above the search threshold).
    # Keep compounds belonging to only one group, even if their score is below the selection threshold.
    score_sel <- 1.1*score_min
    out.clean <- NULL
    for (cmpd in sort(unique(out$Name))) {
        dfident <- out[ out$Name == cmpd, , drop=F]
        if (nrow(dfident)>1) dfident <- dfident[dfident$Score>score_sel, , drop=F]
        if (nrow(dfident)>0) out.clean <- rbind(out.clean, dfident)
    }
    colnames(out.clean) <-  colnames(out)
    out <- as.data.frame(out.clean)
    class(out) <- append(class(out),"annotclusters")
    out
}


#' clusterMerging
#'
#' merging 2 clustering approaches
#'
#' @param clust1 object obtained using the \code{getClusters} function
#' @param clust2 object obtained using the \code{getClusters} function
#' @param dfident1 obtained using the \code{matchClusters} function (optional)
#' @param dfident2 obtained using the \code{matchClusters} function (optional)
clusterMerging <- function(clust1,clust2,dfident1=NULL,dfident2=NULL)
{
    if (! inherits(clust1, "outclust"))
        stop("'clust1' variable is not a 'outclust' class")
    if (! inherits(clust2, "outclust"))
        stop("'clust2' variable is not a 'outclust' class")
    if (!is.null(dfident1) && ! inherits(dfident1, "annotclusters"))
        stop("'dfident1' variable is not a 'annotclusters' class")
    if (!is.null(dfident2) && ! inherits(dfident2, "annotclusters"))
        stop("'dfident2' variable is not a 'annotclusters' class")

   L1 <- sort(unique(dfident1[,2]))
   L2 <- sort(unique(dfident2[,2]))

   NBCLUST <- 0
   clustertab <- NULL
   clusters <- list()
   dfident <- NULL

   # Merging based on annotations
   if (!is.null(dfident1) && !is.null(dfident2)) {
       # List of compounds common to both methods
       for (cmpd in L1[ L1 %in% L2]) {
           ident1 <- dfident1[ which(dfident1[,2]==cmpd), ,drop=F]
           CM1 <- clust1$clustertab[ clust1$clustertab[,2] %in% ident1$Cluster, ,drop=F]
           CM2 <- clust2$clustertab[ clust2$clustertab[,1] %in% CM1[,1], ,drop=F]
           ident2 <- dfident2[ dfident2$Cluster %in% CM2[,2], ,drop=F]
           NBCLUST <- NBCLUST + 1
           if (nrow(ident2)) {
               ctab <- CM2
           } else {
               ctab <- CM1
           }
           clust1$clustertab <- clust1$clustertab[ ! clust1$clustertab[,1] %in% ctab[,1], ,drop=F]
           clust2$clustertab <- clust2$clustertab[ ! clust2$clustertab[,1] %in% ctab[,1], ,drop=F]
       
           dfident1 <- dfident1[ which(dfident1[,2]!=cmpd), ,drop=F]
           dfident2 <- dfident2[ which(dfident2[,2]!=cmpd), ,drop=F]
           #dfident2 <- dfident2[ dfident2$Cluster %in% CM2[,2], ]
       
           idclust <- paste0('C',NBCLUST)
           ctab[,2] <- idclust
           clustertab <- rbind(clustertab, ctab)
           clusters[[idclust]] <- as.numeric(ctab[,3])
           dfident <- rbind(dfident, c(cmpd,idclust))
       }
       
       # List of compounds found only for method 1
       for (cmpd in L1[ ! L1 %in% L2]) {
           ident1 <- dfident1[ which(dfident1[,2]==cmpd), ,drop=F]
           CM1 <- clust1$clustertab[ clust1$clustertab[,2] %in% ident1$Cluster, ,drop=F]
           CM2 <- clust2$clustertab[ clust2$clustertab[,1] %in% CM1[,1], ,drop=F]
           NBCLUST <- NBCLUST + 1
           ctab <- CM1
           clust1$clustertab <- clust1$clustertab[ ! clust1$clustertab[,1] %in% ctab[,1], ,drop=F]
           clust2$clustertab <- clust2$clustertab[ ! clust2$clustertab[,1] %in% ctab[,1], ,drop=F]
       
           dfident1 <- dfident1[ which(dfident1[,2]!=cmpd), ,drop=F]
           dfident2 <- dfident2[ which(dfident2[,2]!=cmpd), ,drop=F]
           #dfident2 <- dfident2[ dfident2$Cluster %in% CM2[,2], ]
       
           idclust <- paste0('C',NBCLUST)
           ctab[,2] <- idclust
           clustertab <- rbind(clustertab, ctab)
           clusters[[idclust]] <- as.numeric(ctab[,3])
           dfident <- rbind(dfident, c(cmpd,idclust))
       }

       ## List of compounds found only for method 2
       for (cmpd in L2[ ! L2 %in% L1]) {
           ident2 <- dfident2[ which(dfident2[,2]==cmpd), ,drop=F]
           CM2 <- clust2$clustertab[ clust2$clustertab[,2] %in% ident2$Cluster, ,drop=F]
           CM1 <- clust1$clustertab[ clust1$clustertab[,1] %in% CM2[,1], ,drop=F]
           NBCLUST <- NBCLUST + 1
           ctab <- CM2
           clust1$clustertab <- clust1$clustertab[ ! clust1$clustertab[,1] %in% ctab[,1], ,drop=F]
           clust2$clustertab <- clust2$clustertab[ ! clust2$clustertab[,1] %in% ctab[,1], ,drop=F]
       
           dfident1 <- dfident1[ which(dfident1[,2]!=cmpd), ,drop=F]
           dfident2 <- dfident2[ which(dfident2[,2]!=cmpd), ,drop=F]
           #dfident2 <- dfident2[ dfident2$Cluster %in% CM2[,2], ]
       
           idclust <- paste0('C',NBCLUST)
           ctab[,2] <- idclust
           clustertab <- rbind(clustertab, ctab)
           clusters[[idclust]] <- as.numeric(ctab[,3])
           dfident <- rbind(dfident, c(cmpd,idclust))
       }
       class(dfident) <- append(class(dfident),"annotclusters")
   }

   ## List of remaining clusters for method 1
   for( CL in unique(sort(clust1$clustertab[,2])) ) {
       CM1 <- clust1$clustertab[clust1$clustertab[,2]==CL, ,drop=F]
       CM2 <- clust2$clustertab[ clust2$clustertab[,1] %in% CM1[,1], ,drop=F]
       NBCLUST <- NBCLUST + 1
       idclust <- paste0('C',NBCLUST)
       CM1[,2] <- idclust
       clustertab <- rbind(clustertab, CM1)
       clusters[[idclust]] <- as.numeric(CM1[,3])
       if (nrow(CM2))
          clust2$clustertab <- clust2$clustertab[ ! clust2$clustertab[,1] %in% CM2[,1], ,drop=F]
   }
   
   ## List of remaining clusters for method 2
   for( CL in unique(sort(clust2$clustertab[,2])) ) {
       CM2 <- clust2$clustertab[clust2$clustertab[,2]==CL, ,drop=F]
       if (c('matrix') %in% class(CM2)) {
          NBCLUST <- NBCLUST + 1
          idclust <- paste0('C',NBCLUST)
          CM2[,2] <- idclust
          clustertab <- rbind(clustertab, CM2)
          clusters[[idclust]] <- as.numeric(CM2[,3])
       }
   }
   class(clusters) <- append(class(clusters), "clusters")

   # Annotated cluster table
   annottab <- NULL
   if (!is.null(dfident)) {
       annottab <- clustertab
       for (i in 1:nrow(dfident))
          annottab[which(annottab[, 2] == dfident[i, 2]), 2] <- dfident[i,1]
   }

   list(clusters=clusters, clustertab=clustertab, dfident=dfident, annottab=annottab, params=list (  method='merging' ))
}
