# ID StatTools.R
# Copyright (C) 2017 INRA
# Authors: D. Jacob
#
trim<-function(x) gsub("^\\s+|\\s+$", "", x)

colors <- c("red", "cornflowerblue", "darkgreen", "blueviolet", "orange", 
             "magenta", "darkred", "coral", "mediumvioletred", "yellow4",
             "seagreen2", "lightskyblue", "darkcyan", "yellowgreen", "limegreen",
             "wheat2", "yellow4",  "violetred1", "darkorange", "cyan4")

matrix <- NULL
factors <- NULL
association <- NULL

options(warn = -1)
options(width=1024)

replace_zero <- function(data)
{
    nbzero <- length(data[data==0])
    minval <- min(data[data>0])
    if (nbzero>0) data[data==0]<- abs( rnorm( nbzero, 0.1*minval, 0.01*minval ) )
    data
}

########################################################################################
## Scaling                                                                            ##
########################################################################################
## Aim: Data scaling
##
## Settings:
##      - matrix: a matrix with the variables in columns and the sample in rows.
##          The first column contains the sample names
##      - the "log" value: According to the data, people can use the logarithmic transformation.
##          If log=0: no logaritmic transformation
##          If log=2: data are transformed with log2
##          If log=10: data are transformed with log10
##      - Scaling methods: Several methods can be applied to a data set. In this application, we put forward:
##          Centering
##          Zscore
##          Pareto
##          Vast
##          Range
##          Level
##          L1: (Euclidean norm)
##          L2: (Euclidean norm)
##          NULL
##
########################################################################################
get_Scaling <- function(matrix, log, methods=c("Centering","Zscore","Pareto","Vast","Range","Level","L1","L2","NULL"))
{
## LOG Transformation
   if (log==2) {
       matrix <- replace_zero(matrix)
       matrix<-log2((matrix))
       cat('transformation to log2\n')
   }
   if (log==10) {
       matrix <- replace_zero(matrix)
       matrix<-log10(matrix)
       cat('transformation to log10\n')
   }
   if (log==11) {
       matrix<-log10(matrix + 1)
       cat('Transformation to log10 after incrementing (+1)\n')
   }

## Column-wise
   if ("Centering" %in% methods)   { matrix<-scale(matrix,center=T,scale=F) }
   if ("Zscore" %in% methods)      { matrix<-scale(matrix,center=T,scale=T) }
   if ("Pareto" %in% methods)      { matrix<-apply(t(matrix), 1,function(x){(x-mean(x))/sqrt(sd(x))}) }
   if ("Vast" %in% methods)        { matrix<-apply(t(matrix), 1,function(x){(x-mean(x))*mean(x)/(sd(x)*sd(x))}) }
   if ("Range" %in% methods)       { matrix<-apply(t(matrix), 1,function(x){(x-mean(x))/(max(x)-min(x))}) }
   if ("Level" %in% methods)       { matrix<-apply(t(matrix), 1,function(x){(x-mean(x))/median(x)}) }

## Row-wise
   if ("L1" %in% methods)          { matrix<-scale(matrix,center=T,scale=F)
                                     matrix<-t(apply(matrix, 1,function(x){x/sum(x)})) }
   if ("L2" %in% methods)          { matrix<-scale(matrix,center=T,scale=F)
                                     matrix<-t(apply(matrix, 1,function(x){x/sqrt(sum(x*x))}))  }

   return(matrix)
}


#====================================================================#
# Clustering - Find the optimal Cutting Tree Value
#====================================================================#

get_Clusters <- function(data, method='hca', ... )
{
    out <- NULL
    if (method=='hca')  out <- get_Clusters_hca(data, ... )
    if (method=='corr') out <- get_Clusters_corr(data, ... )
    out
}

get_Clusters_hca <- function(data, vcutusr=0, scalemeth='Zscore', log=0, distmeth='euclidean', hcmeth='complete', maxcsize=40, mincsize=2, bucchar='B')
{
   params <- list (
     method='hca',

   # Scaling params
     scalemeth = scalemeth, 
     log = log,

   # Bucket variables
     BUCCHAR = bucchar,
   
   # Clustering
     DISTMETH = distmeth,
     HCMETH = hcmeth,
     MAXCSIZE = maxcsize,
     MINCSIZE = mincsize,

   # the range of 'Vcut' for getting some statistics
     Vcut_min=0.1,
     Vcut_max=0.9,
     Vstep=0.005,
     vcutusr=vcutusr

   )

   # Get the association table: Variables <=> CLusterID,
   # by applying the cutting tree threshold 'vcut' on the Hclust 'hc'
   get_VarSet <- function(hc, vcut)
   {
        hcut <- vcut*( max(hc$height) + min(hc$height) )
        CT <- cutree(hc.c,h=hcut)
        i <- 1
        V <- NULL
        for (k in 1:length(unique(CT))) {
            if (length(CT[CT==k])<params$MINCSIZE) next
            CTK <- CT[CT==k]
            CTK[] <- i
            V <- c(V, CTK)
            i <- i + 1
        }
        VarSet <- cbind( names(V), simplify2array(lapply( as.vector(V), function(x) { paste('C',x,sep=''); })) )
        return(VarSet)
   }

   matrix <- get_Scaling(data,log=params$log, methods=params$scalemeth)

   VC <- 1
   CC <- 2
   Vcut_min=params$Vcut_min
   Vcut_max=params$Vcut_max
   Vstep=params$Vstep

# Variables Classification
   x <- as.matrix(matrix)

# Compute the distance matrix (see ?dist) and the Hierarchical clustering (see ?hclust)
   dist.c <- dist(t(x), method=params$DISTMETH)
   hc.c <- hclust(dist.c, method=params$HCMETH)

# Init
   Vstats <- NULL
   Vcut <- seq(from=Vcut_min, to=Vcut_max, by=Vstep)

# Loop on the range
   for( k in 1:length(Vcut) ) {

     # Init
     Cluster_Count <- 0
     Cluster_VarCount <- 0
     Cluster_MatchVarCount <- 0
     Cluster_MaxSize <- 0
     Cluster_Matching <- 0
     VarSet <- get_VarSet(hc.c, Vcut[k])
     Clusters <- as.character(unique(VarSet[,2]))

     # For each cluster ...
     for (i in 1:length(Clusters)) {
        CL <- Clusters[i]
        Cluster_Count <- Cluster_Count + 1

        # get the peaklist
        ppmlst  <- VarSet[VarSet[,CC] == CL, VC]
        Cluster_VarCount <- Cluster_VarCount + length(ppmlst)
        if (Cluster_MaxSize < length(ppmlst) ) Cluster_MaxSize <- length(ppmlst)
     }
     Vstats <- rbind( Vstats, c(Vcut[k], Cluster_Count, Cluster_VarCount, Cluster_MaxSize ) )

   }
   colnames(Vstats) <- c("Vcut", "Nb Clusters", "Nb Vars", "MaxSize")

# Compute the optimal value for Vcut
# -- find the set of Vcut values such as the cluster count is max
   Vsub <- Vstats[ Vstats[,4]<params$MAXCSIZE,  ]
   V <- Vsub[Vsub[,2] == max(Vsub[,2]),1]
# -- Keep the Vcut value such as the variable count is max
   VcutOpt <- V[which(Vstats[Vstats[,1] %in% V,3]==max(Vstats[Vstats[,1] %in% V,3]))][1]
   Vcut <- ifelse(params$vcutusr==0, VcutOpt, params$vcutusr)
   indx <- round(Vcut/Vstep) - round(Vcut_min/Vstep) + 1

# Get the association table: Variables <=> CLusterID
   association <- get_VarSet(hc.c, Vcut)
   association <- cbind( association, gsub(params$BUCCHAR, "", gsub("_", ".", association[,VC])) )

   lclust <- unique(sort(association[,2]))
   clusters <- list()
   for(i in 1:length(lclust))  {
       CL <- lclust[i]
       clusters[[CL]] <-  as.numeric(association[association[,2]==CL, ][,3])
   }

   cat("#-- Clustering --\n")
   cat('#  Distance Method:',params$DISTMETH,"\n")
   cat('#  Agglomeration Method:',params$HCMETH,"\n")
   cat('#  Cutting Tree threshold:',Vcut,"\n")
   cat('#  Nb Clusters:',length(lclust),"\n")
   cat("#\n")

   return( list(vstats=Vstats, clusters=clusters, clustertab=association, params=params, vcrit=Vcut, indxopt=indx) )

}

get_Clusters_corr <- function(data, scalemeth='Zscore', log=0, cmeth='pearson', cval=0, dC=0.01, maxcsize=40, mincsize=2, bucchar='B', stats_comp=T, ncpu=1)
{
   params <- list (
     method='corr',

   # Scaling params
     scalemeth = scalemeth, 
     log = log,

   # Bucket variables
     BUCCHAR = bucchar,

   # Clustering
     CMETH = cmeth,
     CVAL = cval,
     dC = dC,
     MINCSIZE = mincsize,
     MAXCSIZE = maxcsize,

   # get some statistics
     VSTATS_CAL=stats_comp,
     NCPU=ncpu,

   # the range of 'Cval' for getting some statistics
     CVAL_MIN=0.9,
     CVAL_MAX=0.999,
     CVAL_STEP=0.001

   )

   matrix <- get_Scaling(data,log=params$log, methods=params$scalemeth)

#---- CLUSTERING -----

   require(igraph)

 # Method in ("pearson", "kendall", "spearman")
   cor_mat <- cor(matrix,method=params$CMETH)
   cor_mat[ lower.tri(cor_mat, diag=TRUE) ] <- 0

   vstats <- NULL
   indx <- 0
   if (params$VSTATS_CAL || params$CVAL==0) {
       cvals <- seq(from=params$CVAL_MIN, to=params$CVAL_MAX, by=params$CVAL_STEP)
       vstats <- estime_cval(cor_mat, cvals, params$NCPU)
       if (params$CVAL==0) {
           sub1 <- which( vstats[,4]<params$MAXCSIZE )
           sub2 <- which ( vstats[ sub1, 2 ]==max(vstats[ sub1, 2 ] ) ) + sub1[1] - 1
           indx <- sub2[ which(vstats[sub2,4]==min(vstats[sub2,4]))[1] ]
           params$CVAL <- vstats[indx, 1]
       } else {
           indx <- round(params$CVAL/params$CVAL_STEP) - round(params$CVAL_MIN/params$CVAL_STEP) + 1
       }
   }

   M <- NULL
   M <- cbind(M,colnames(matrix))
   NVARS <- length(colnames(matrix))

   cvalset <- c(params$CVAL-params$dC,params$CVAL,params$CVAL+params$dC)
   NBCLUST <- vector(mode="integer", length=length(cvalset))
   for (k in 1:length(cvalset)) {
      cval <- cvalset[k]
      cor_mat_cval <- cor_mat
      cor_mat_cval[ cor_mat_cval < cval] <- 0
      M<-cbind(M,c(1:NVARS)*0)
      graph <- graph.adjacency(cor_mat_cval>cval, weighted=TRUE, mode="upper")
      E(graph)$weight<-t(cor_mat_cval)[t(cor_mat_cval)>cval]
      V(graph)$label<- V(graph)$name
      cliques <- sapply(decompose.graph(graph), vcount)
      ordcli <- order(cliques,decreasing = T)
      g<-0
      for (i in 1:length(ordcli)) {
          ind<-ordcli[i]
          if (cliques[ind]>=params$MINCSIZE) {
             g <- g + 1
             subg <- decompose.graph(graph)[[ind]]
             VARS<-V(subg)$name
             M[M[,1] %in% VARS,k+1] <- sprintf("C%d",g)
          }
      }
      NBCLUST[k] <- g
   }

   MC <- NULL
   for (i in 1:NVARS) {
      if (is.na(sum(as.numeric(M[i,-1])))) MC <- rbind(MC,M[i,])
   }
   varnames <- as.vector(lapply(cvalset, function(x) sprintf('V_%s',x,2)))
   varnames <- c('ppm',varnames)
   colnames(MC) <- varnames

   # Synthesis
   for (s in c(1,2)) {
      for (k in 1:NBCLUST[s]) {
         V <- MC[MC[,1+s]==sprintf("C%d",k),c(1,1+s,2+s)]
         i<-1
         while(i<dim(V)[1]) {
            if (V[i,3] == "0") {
               n<-1
               while(i<dim(V)[1] && V[i+1,2]==V[i,2] && V[i+1,3] == "0") { i<-i+1; n<-n+1; }
               # Rule 1
               if (n==1) {
                  if (i>1 && V[i-1,2]==V[i,2]) V[i,3] <- V[i-1,3]
                  else if (i<dim(V)[1] && V[i+1,2]==V[i,2]) V[i,3] <- V[i+1,3];
               }
               # Rule 2
               if (n>1) {
                  NBCLUST[1+s] <- NBCLUST[1+s] + 1
                  for (j in 1:n) V[i-n+j,3] <- sprintf("C%d",NBCLUST[1+s])
               }
            }
            i<-i+1
         }
         MC[MC[,1] %in% V[,1],2+s] <- V[,3]
      }
   }
   association <- MC[,c(1,4)]
   association <- association[association[,2] != "0",]
   association <- cbind( association, gsub(params$BUCCHAR, "", gsub("_", ".", association[,1])) )
   colnames(association) <- c("VAR","CLID", "PPM")

   lclust <- unique(sort(association[,2]))
   clusters <- list()
   for(i in 1:length(lclust))  {
       CL <- lclust[i]
       clusters[[CL]] <-  as.numeric(association[association[,2]==CL, ][,3])
   }

   cat("#-- Clustering --\n")
   cat('#  Correlation Method:',params$CMETH,"\n")
   cat('#  Correlation Threshold :',params$CVAL,"\n")
   cat('#  Correlation Tolerance:',params$dC,"\n")
   cat('#  Nb Clusters:',length(lclust),"\n")
   cat("#\n")

   return( list( vstats=vstats, clusters=clusters, clustertab=association, params=params, vcrit=params$CVAL, indxopt=indx ) )
}

estime_cval <- function(cor_mat, cvals, ncpu)
{
   require(igraph)
   require(doParallel)

   compute_crit <- function(cor_mat, cval) {
      cor_mati <- cor_mat
      cor_mati[ cor_mati < cval] <- 0
   
      graph <- graph.adjacency(cor_mati>cval, weighted=TRUE, mode="upper")
      E(graph)$weight<-t(cor_mati)[t(cor_mati)>cval]
      V(graph)$label<- V(graph)$name
      cliques <- sapply(decompose.graph(graph), vcount)
      ordcli <- order(cliques,decreasing = T)
      nb_clusters<-0
      nb_clusters_2<-0
      nb_buckets<-0
      for (i in 1:length(ordcli)) {
         if (cliques[ordcli[i]]>=2) {
           nb_clusters <- nb_clusters + 1
           nb_buckets <- nb_buckets + cliques[ordcli[i]]
         }
         if (cliques[ordcli[i]]==2)
           nb_clusters_2 <- nb_clusters_2 + 1
      }
      size_max <- cliques[ordcli[1]]
      CRIT<-size_max/nb_clusters
      c(cval,nb_clusters,nb_buckets,size_max,nb_clusters_2,-20*log10(CRIT))
   }
   
   cl <- makeCluster(ncpu)
   registerDoParallel(cl)

   vstats <- foreach(i=1:length(cvals), .combine=rbind, .packages = c("igraph", "doParallel")) %dopar% {  compute_crit(cor_mat, cvals[i]); }

   stopCluster(cl)

   colnames(vstats) <- c("Cval","Nb Clusters","Nb Vars","Max Size","Nb Clusters 2","Criterion")
   vstats
}


get_Merged_dataset <- function(data, clustObj, onlycluster=FALSE)
{
   association <- clustObj$clustertab
   matrix <- data
   if (class(association)=="matrix" && dim(association)[1]>1) {
      Clusters <- unique(association[,2])
      M <- simplify2array(lapply( 1:length(Clusters), function(x) { apply( matrix[, colnames(matrix) %in% association[association[,2] == Clusters[x],1] ], 1, mean ) } ))
      colnames(M) <- Clusters
      if (!onlycluster)
         M <- cbind(M, matrix[, !(colnames(matrix) %in% association[,1])])
      matrix <- M
   }
   matrix
}

plot.Criterion <- function(clustObj)
{
    Vstats <- clustObj$vstats
    Vcrit <- clustObj$vcrit
    n <- clustObj$indxopt
    params <- clustObj$params
    MaxSize <- params$MAXCSIZE
    method <- params$method

   # PNG image - Vstats plots
     legNames <- c( colnames(Vstats)[3], colnames(Vstats)[2], colnames(Vstats)[4])
     par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
     xlim <- c(min(Vstats[,1]), max(Vstats[,1]))
     plot(Vstats[,1],Vstats[,3], xlim=xlim, ylim=c(0,1.2*max(Vstats[,3])), type="l", col=colors[1],
         xlab="Critere", ylab="Number of variables / Clust Max Size",
         main=sprintf("Critere = %4.3f,  Nb Clust = %d,  Nb Vars = %d,  Clust Max Size = %d",Vcrit, Vstats[n,2], Vstats[n,3], Vstats[n,4]))
     lines(Vstats[,1],Vstats[,4], col=colors[3])
     abline(v=Vcrit, col="red")
     abline(h=MaxSize, col="green")
     if( method=='hca' ) {
          slines <- seq(from=params$Vcut_min, to=params$Vcut_max, by=0.025)
     } else { 
          slines <- seq(from=params$CVAL_MIN, to=params$CVAL_MAX, by=0.0025)
     }
     abline(v = slines, col = "lightgray", lty = 3)
     legend("topright", legNames, col=colors[c(1:length(legNames))], lwd=2, lty=1, cex=0.6, horiz=FALSE)
     par(new=TRUE)
     plot(Vstats[,1],Vstats[,2], xlim=xlim, type="l", col=colors[2], xaxt="n", yaxt="n", xlab="",ylab="" )
     axis(side=4, at = pretty(range(Vstats[,2])))
     mtext("Number of Clusters", side=4, line=3)
}

plot.Clusters <- function(data, clustObj)
{
     CLUST <- unique(clustObj$clustertab[,2])
     X <- get_Merged_dataset(data, outclust, onlycluster=T)
     X[X<=0]<-0.00001*max(X)
     cent <- log10(t(X))

     Rank <- simplify2array(lapply(c(1:length(CLUST)), function(x) { round(mean(cent[x, ], na.rm=T))  }))
     cols <- c('red', 'orange', 'darkgreen', 'blue', 'purple')
     Rank <- Rank  - min(Rank) + 1
     nbrank <- max(Rank) - min(Rank) + 1

     boxplot(t(cent),cex.axis=0.5, main="Boxplot by clusters (log10 transformed)",
          ylab='\n\n\n\n\nClusters',xlab='Intensity (log10)', outline=F, horizontal=T, border=cols[Rank], las=2, cex.axis=0.5)
     abline(v=c(1:length(CLUST)), col="gray95")
     par(new=TRUE)
     boxplot(t(cent),cex.axis=0.5, main="Boxplot by clusters (log10 transformed)",
          ylab='\n\n\n\n\nClusters',xlab='Intensity (log10)', outline=F, horizontal=T, border=cols[Rank], las=2, cex.axis=0.5)
}

##########################################################################
## Fonction plot.loadings                                               ##
##########################################################################
## Loadings plot : Input : data [N,PC] N variables (rows), PC (columns)
## 
## 
## use with associations file:
## 
## Note : the size of the ellipses is related to the level of confidence chosen p: 
## p = 1 - exp(-cellipse^2 / 2) => cellipse <- sqrt(qchisq(p, 2))
## With a confidence level at 66.66% (2/3), we have cellipse = 1.5
## See http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
##     http://www.earth-time.org/projects/upb/public_docs/ErrorEllipses.pdf
##########################################################################

plot.Loadings <- function (data,pc1,pc2, associations=NULL,
         main="Loadings", xlimu=c(min(data[,pc1]),max(data[,pc1])), ylimu=c(min(data[,pc2]),max(data[,pc2])), cexlabel=1)
{
   V <- t(data[,c(pc1,pc2)])
   colnames(V) <- rownames(data)
   plot(t(V),pch=20, adj=1, lwd=1, main=main,
            xlab= colnames(data)[pc1], ylab= colnames(data)[pc2], 
            xlim=xlimu, ylim=ylimu, col='red')
   abline(h=0,v=0)
   plot.ellipse( data[, pc1], data[, pc2], center=c(0,0), level=0.666, col="red", lty=3, lwd=0.1, type="l")

   # Use Associations (assoc1=T)
   if (class(associations)=="matrix" && dim(associations)[1]>1) {
          points(t(V), pch=19, col="grey")
          text( adj=0, t(V)[,1], t(V)[,2], col="lightgrey", colnames(V), cex=cexlabel )
          colnames(V) <- as.vector(sapply(rownames(data),function(x) { ifelse(sum( associations[,1] %in% x), associations[,2][associations[,1] %in%  x ], x) }))
          P <- V[,colnames(V) %in% associations[,2] ]
          cols <- "red"
          Clusters <- unique(associations[order(associations[,1],decreasing=T),2])
          if (length(Clusters)<length(colnames(P))) {
             for (i in 1:length(Clusters)) {
                XY <- t(P[,colnames(P)==Clusters[i]])
                M<-c(mean(XY[,1]),mean(XY[,2]))
                if (dim(XY)[1]>1) {
                   if (dim(XY)[1]>2)
  				      plot.ellipse( XY[,1], XY[,2], center=M, level= 0.8646647, col=colors[i %% length(colors)], lty=3, lwd=0.1, type="l")
                   else
                      for (j in 1:dim(XY)[1]) lines(rbind( XY[j, ] , M), col="black")
                }
                points(XY, pch=19, col=colors[i %% length(colors)])
                text(adj=0, M[1], M[2], Clusters[i], col="black", cex=cexlabel)
             }
             for (i in 1:length(Clusters)) {
                XY <- t(P[,colnames(P)==Clusters[i]])
                M<-c(mean(XY[,1]),mean(XY[,2]))
                text(adj=0, M[1], M[2], Clusters[i], col="black", cex=cexlabel)
             }
          }
          else {
             col_ids <- c(1:length(Clusters))
             cols <- NULL; for (i in 1:length(colnames(P))) cols <- c(cols,  col_ids[Clusters == colnames(P)[i]])
             points(t(P), pch=19, col="red")
             text( adj=0, t(P)[,1], t(P)[,2], col=cols, colnames(P), cex=cexlabel )
          }
   }
   else {
      text( adj=0, t(V)[,1], t(V)[,2], col="cornflowerblue", colnames(V), cex=cexlabel )
   }
}

##########################################################################
## Fonction plot.loadings                                               ##
##########################################################################
## Compute confidence ellipses.
##  x and y variables for drawing.
##  level confidence level used to construct the ellipses. By default, 0.95.
##  npoint number of points used to draw the ellipses.
##  bary logical value. If TRUE, the coordinates of the ellipse around the barycentre of individuals are calculated.
## See https://github.com/kassambara/ggpubr/blob/master/R/stat_conf_ellipse.R
##########################################################################

plot.ellipse <- function ( x, y, center=NULL, level = 0.95, npoint = 100, bary = FALSE, ... )
{
  .ellipse <- function(x, scale = c(1, 1), centre = c(0, 0), t = 1.5, which = c(1,2), npoints = 100) {
    names <- c("x", "y")
    if (is.matrix(x)) {
      xind <- which[1]
      yind <- which[2]
      r <- x[xind, yind]
      if (missing(scale)) {
        scale <- sqrt(c(x[xind, xind], x[yind, yind]))
        if (scale[1] > 0)
          r <- r/scale[1]
        if (scale[2] > 0)
          r <- r/scale[2]
      }
      if (!is.null(dimnames(x)[[1]]))
        names <- dimnames(x)[[1]][c(xind, yind)]
    }
    else r <- x
    r <- min(max(r, -1), 1)
    d <- acos(r)
    a <- seq(0, 2 * pi, len = npoints)
    matrix(c(t*scale[1]*cos(a + d/2) + centre[1], t*scale[2]*cos(a - d/2) + centre[2]), npoints, 2, dimnames = list(NULL, names))
  }


  if (is.null(center)) center <- c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))
  tab <- data.frame(x = x, y = y)
  mat.cov <- stats::cov(tab)
  t = sqrt(stats::qchisq(level, 2))
  if (bary)
  mat.cov = mat.cov/nrow(tab)
  res <- .ellipse(mat.cov, centre = center, t=t, npoints = npoint)
  lines(res, adj=1, main="", xlab="", ylab="", ... )
}
