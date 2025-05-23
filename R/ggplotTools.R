# ID ggplotTools.R
# Copyright (C) 2019 INRA
# Authors: D. Jacob
#

#-------------------------------------------------------------------------------------------------------------------------
# ggplot2 routines
#-------------------------------------------------------------------------------------------------------------------------

##########################################################################
## See https://github.com/hrbrmstr/ggalt/blob/master/R/geom_encircle.r
## author: Ben Bolker
##########################################################################

draw_key_hack <- function(data, params, size) {
  data$fill <- scales::alpha(data$fill, data$alpha)
  data$alpha <- 1

  grid::grobTree(
    if (!is.na(data$fill)) grid::rectGrob(gp = grid::gpar(col = NA, fill = data$fill)),
    ggplot2::draw_key_path(data, params)
  )
}

GeomEncircle <- ggplot2::ggproto("GeomEncircle", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(colour = "black",
                    fill = NA, ## ???
                    alpha = 1,
                    linetype=1,
                    size=1,
                    s_shape=0.5,  ## corresponds to default shape in xspline of -0.5
                    s_open=FALSE,
                    expand=0.05,
                    spread=0.1),
  draw_key = draw_key_hack, ## ???

  draw_group = function(data, panel_scales, coord) {
    ## browser()
    coords <- coord$transform(data, panel_scales)
    first_row <- coords[1, , drop = FALSE]
    rownames(first_row) <- NULL ## prevent warning later

    m <- lapply(coords[,c("x","y")],mean,na.rm=TRUE)
    ch <- grDevices::chull(coords[c("x","y")])

    mkcoords <- function(x,y) {
        data.frame(x,y,first_row[!names(first_row) %in% c("x","y")])
    }

    coords <- coords[ch,]
    ## FIXME: using grid:: a lot. importFrom instead?

    ## convert from lengths to physical units, for computing *directions*
    cc <- function(x,dir="x")
        grid::convertUnit(grid::unit(x,"native"),"mm",typeFrom="dimension",
                          axisFrom=dir,valueOnly=TRUE)

    ## convert back to native (e.g. native + snpc offset)
    cc_inv <- function(x,dir="x")
        grid::convertUnit(x,"native",typeFrom="location",
                          axisFrom=dir,valueOnly=TRUE)

    cc_comb <- function(x1,x2,dir="x")
        cc_inv(unit(x1,"native")+unit(x2,"snpc"),dir=dir)

    ## find normalized vector: d1 and d2 have $x, $y elements
    normFun <- function(d1,d2) {
        dx <- cc(d1$x-d2$x)
        dy <- cc(d1$y-d2$y)
        r <- sqrt(dx*dx+dy*dy)
        list(x=dx/r,y=dy/r)
    }

    if (nrow(coords)==1) {
        ## only one point: make a diamond by spreading points vertically
        ## and horizontally
        coords <- with(coords,
                       mkcoords(
                           c(x,x+spread,x,x-spread),
                           c(y+spread,y,y-spread,y)))
    } else if (nrow(coords)==2) {
        ## only two points: make a diamond by spreading points perpendicularly
        rot <- matrix(c(0,1,-1,0),2)
        dd <- c(rot %*% unlist(normFun(coords[1,],coords[2,])))*
                     coords$spread
        coords <- with(coords, {
            ## figure out rotated values, then convert *back* to native units
            ## already in scaled units, so ignore?
            x <- c(x[1],
                   m$x+dd[1], ## cc_comb(m$x,dd[1]),
                   x[2],
                   m$x-dd[1]) ## cc_comb(m$x,-dd[1]))
            y <- c(y[1],
                   m$y+dd[2], ## cc_comb(m$y,dd[2],"y"),
                   y[2],
                   m$y-dd[2]) ## cc_comb(m$y,-dd[2],"y"))
            mkcoords(x,y)
        })
    }

    disp <- normFun(coords,m)

    ## browser()

    gp <- grid::get.gpar()
    pars1 <- c("colour","linetype","alpha","fill","size")
    pars2 <- c("col","lty","alpha","fill","lwd")
    gp[pars2] <- first_row[pars1]
    grid::xsplineGrob(
        with(coords,unit(x,"npc")+disp$x*unit(expand,"snpc")),
        with(coords,unit(y,"npc")+disp$y*unit(expand,"snpc")),
        ## coords$x,
        ## coords$y,
      shape = coords$s_shape-1,  ## kluge!
      open = first_row$s_open,
      gp = gp)
  }
)

##########################################################################
## Fonction geom_encircle
##########################################################################
geom_encircle <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomEncircle, mapping = mapping,  data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' ggplotCriterion
#'
#' Plots the curves that show the number of clusters, the number of clustered buckets and the 
#' size of biggest cluster  versus the criterion, namely the correlation threshold for the 'corr' 
#' method, the cutting value for the 'hca' method.
#'
#' @param clustObj a list generated by the \code{getClusters} function
#' @param reverse indicates if the x axis need to be reversed
ggplotCriterion <- function(clustObj, reverse=FALSE)
{
    Vstats <- clustObj$vstats
    Vcrit <- clustObj$vcrit
    n <- clustObj$indxopt
    params <- clustObj$params
    MaxSize <- params$MAXCSIZE
    method <- params$method
    associations <- clustObj$clustertab
    lclust <- unique(sort(associations[,2]))

    dC <- 0
    if (method=="corr") dC <- params$dC

    dfstat <- data.frame(Vcut=Vstats[,1], nbclust=Vstats[,2], nbvars=Vstats[,3], Maxsize=Vstats[,4])

    ratioVC <- round(max(dfstat$nbvars)/max(dfstat$nbclust),2)

    GTITLE <- sprintf("%s method, Critere = %4.3f,  Nb Clust = %d",toupper(method),Vcrit, length(lclust))
    xlab <- "Critere"
    ylab <- paste0("Nb of variables / Clust Max Size / Nb of Cluster x ",ratioVC)

    Vcut <- nbclust <- nbvars <- Maxsize <- NULL
    g <- ggplot2::ggplot(dfstat, ggplot2::aes(x = Vcut)) + ggplot2::ggtitle(GTITLE) +
         ggplot2::geom_line(ggplot2::aes(y = nbclust*ratioVC, colour="nbclust")) +
         ggplot2::geom_line(ggplot2::aes(y = nbvars, colour="nbvars")) +
         ggplot2::geom_line(ggplot2::aes(y = Maxsize, colour="Maxsize")) +
         ggplot2::geom_vline(xintercept=Vcrit-dC, color="grey", linetype="dashed", , size=0.2) +
         ggplot2::geom_vline(xintercept=Vcrit+dC, color="grey", linetype="dashed", , size=0.2) +
         ggplot2::geom_vline(xintercept=Vcrit, color="red", linetype="dashed", , size=0.2) +
         ggplot2::annotate("rect", xmin=Vcrit-dC, xmax=Vcrit+dC, ymin=0, ymax=Inf, alpha=0.2, fill="grey") +
         ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~.*1/ratioVC, name = "Nb Clust")) +
         ggplot2::scale_colour_manual("", breaks = c("nbvars", "nbclust", "Maxsize"), values = c( "darkgreen", "blue", "red")) +
         ggplot2::labs(y = ylab, x = xlab) +
         ggplot2::theme_light() + ggplot2::theme(legend.position="bottom")
    if (max(Vstats[,2])>MaxSize) g <- g + ggplot2::geom_hline(yintercept=MaxSize, color="green", linetype="dashed", , size=0.2)
    if (reverse) g <- g + ggplot2::scale_x_reverse()
    g
}


#' ggplotClusters
#'
#' Plots the boxplot of all clusters allowing to have an insight on the clusters distribution.
#' Plot based on ggplot2
#'
#' @param data the matrix including the integrations of the areas defined by the buckets (columns)
#' on each spectrum (rows)
#' @param clustObj a list generated by the \code{getClusters} function
ggplotClusters <- function(data, clustObj)
{
     CLUST <- unique(clustObj$clustertab[,2])
     X <- getMergedDataset(data, clustObj, onlycluster=T)
     X[X<=0]<-0.00001*max(X)
     cent <- log10(t(X))

     Rank <- simplify2array(lapply(c(1:length(CLUST)), function(x) { round(mean(cent[x, ], na.rm=T))  }))
     cols <- c('red', 'orange', 'darkgreen', 'blue', 'purple')
     Rank <- Rank  - min(Rank) + 1
     nbrank <- max(Rank) - min(Rank) + 1

     Rank <- simplify2array(lapply(c(1:length(CLUST)), function(x) { round(mean(cent[x, ], na.rm=T))  }))
     Rank <- Rank  - min(Rank) + 1

     M <- cbind( Rank, rownames(cent), cent)
     colnames(M)[1:2] <- c("Rank", "Clusters")
     M <- as.data.frame(M, stringsAsFactors=FALSE)
     N <- stats::reshape(M, direction="long", idvar=c("Rank","Clusters"), varying=list(3:dim(M)[2]))

     df <- data.frame( Group=as.factor(N[,1]), Clusters=as.factor(N[,2]), Values=as.numeric(N[, 4]) )

     Clusters <- Values <- Group <- NULL
     method <- clustObj$params$method
     g <- ggplot2::ggplot(df, ggplot2::aes(x=Clusters,y=Values)) + ggplot2::geom_boxplot(ggplot2::aes(fill=Group), outlier.colour = "red", outlier.shape = 1) +
          ggplot2::labs(y = "Intensity (log10)", x = "Clusters") + ggplot2::ggtitle(paste0("Boxplot by clusters (",toupper(method),")")) +
          ggplot2::theme_light() + ggplot2::theme(legend.position="none")
     g
}

# geom_cluster - internal routines
#    Adds a figure layer (polygons or ellipse) grouping the identical identifiers  - based on the CLID field of the dataframe 'data' 
# params
#    g : ggplot2 object
#    data : input dataframe
#    draw.contour: type of contour; possible values are : 'ellipse', 'polygon', 'ellipse2', 'none'
#    min.size : Minimum size of a set of points with the same CLID to draw the corresponding contour
#    npoints : number of points for drawing each contour
#    level : confidence level in case where 'ellipse' is chosen as contour
#    lw: line width, ps: point size, fs: font size
#    draw.points : if TRUE, draw the points
#    draw.labels : if TRUE, draw the labels
#    color.labels : specify a vector of colors for contours (NULL by default)
geom_cluster <- function(g=NULL, data, level=0.8, lw=0.3, ps=0.5, fs=3, min.size=2, npoints=50, 
                         draw.contour="ellipse", draw.points=TRUE, draw.labels=TRUE, color.labels=NULL)
{
    ## Compute confidence ellipses.
    ##  x covariance matrix
    ##  level confidence level used to construct the ellipses. By default, 0.95.
    ##  npoint number of points used to draw the ellipses.
    .ellipse <- function(x, scale = c(1, 1), centre = c(0, 0), level = 0.95, which = c(1,2), npoints = 100)
    {
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
        t = sqrt(stats::qchisq(level, 2))
        r <- min(max(r, -1), 1)
        d <- acos(r)
        a <- seq(0, 2 * pi, len = npoints)
        matrix(c(t*scale[1]*cos(a + d/2) + centre[1], t*scale[2]*cos(a - d/2) + centre[2]), npoints, 2, dimnames = list(NULL, names))
    }
    expand_polygon <- function(X,expand=0.1,tol=1e-3)
    {
       X2 <- X
       nl <- dim(X)[1]
       for( i in 1:nl ) {
           if(i==1)        idx <- c( nl,  i+1)
           if(i==nl)       idx <- c( i-1, 1 )
           if(i<nl && i>1) idx <- c( i-1, i+1 )
           dx <- X[i,1] - 0.5*(X[ idx[1], 1] + X[ idx[2], 1])
           dy <- X[i,2] - 0.5*(X[ idx[1], 2] + X[ idx[2], 2])
           dom <- sqrt((dx)^2 + (dy)^2)
           if ( abs(dx)<tol ) { X2[i, 1:2] <- c(X[i,1], X[i,2] + sign(dy)*expand*dom);  next; }
           if ( abs(dy)<tol ) { X2[i, 1:2] <- c(X[i,1] + sign(dx)*expand*dom, X[i,2] ); next; }
           m=dy/dx; delta <- expand*dom/sqrt(m*m+1)
           X2[i, 1:2] <- c( X[i,1] + sign(dx)*delta, X[i,2] + sign(dy)*sign(m)*m*delta );
       }
       X2
    }

    clustids <- unique(data$CLID)
    clustsize <- sapply(1:length(clustids), function(x) { sum(data$CLID==clustids[x]) })
    dfsub <- data[data$CLID %in% clustids[clustsize>=min.size],]
    subCL <- sort(unique(dfsub$CLID))
    centroids <- stats::aggregate(cbind(pc1,pc2) ~ CLID , dfsub, mean)
    centroids$color=grDevices::rainbow(length(subCL), s=0.85, v=0.7)
    dfsub$color <- sapply( dfsub$CLID, function(x) { centroids[centroids$CLID==x,]$color })

    pc1 <- pc2 <- CLID <- NULL
    for( i in 1:length(subCL) ) {
       dfcl <- dfsub[dfsub$CLID==subCL[i],]
       dfct <- centroids[centroids$CLID==subCL[i],]
       if (sum(subCL[i] %in% clustids[clustsize==2])) {
          if (draw.contour!="none")
              g <- g + ggplot2::geom_path(data=dfcl, color=dfct$color, alpha = 0.2)
       } else {
          if (sum(grep("ellipse", draw.contour))>0) {
              dfel <- data.frame(CLID=subCL[i], .ellipse(stats::cov(dfcl[, c("pc1","pc2")]), 
                                 centre=as.matrix(dfct[,c("pc1","pc2")]), level=level, npoints=npoints),stringsAsFactors=FALSE)
              if (draw.contour=="ellipse")
                 g <- g + ggplot2::geom_polygon(data = dfel, ggplot2::aes(x=pc1, y=pc2), fill=dfct$color, alpha = 0.2, show.legend=FALSE)
              if (draw.contour=="ellipse2")
                 g <- g + ggplot2::geom_path(data=dfel, color=dfct$color, size=lw, alpha = 0.5)
          }
          if (draw.contour=="polygon") {
              dfel <- expand_polygon(dfcl[grDevices::chull(dfcl[, 1:2]), ], expand=0.2)
              g <- g + ggplot2::geom_polygon(data = dfel, ggplot2::aes(x=pc1, y=pc2), fill=dfct$color, alpha = 0.2)
          }
          if (draw.contour=="encircle") # library(ggalt)
              g <- g + geom_encircle(data=dfcl, ggplot2::aes(x=pc1, y=pc2), fill=dfct$color, expand=0.01, s_shape=0.1, spread=0.04, alpha=0.2, colour="white")

       }
       if (draw.points) g <- g + ggplot2::geom_point(ggplot2::aes(x=pc1, y=pc2), data=dfcl, color=dfct$color, size=ps)
   }
   if (draw.labels) {
         thecolors <- centroids$color
         if (! is.null(color.labels)) thecolors <- color.labels
         g <- g + ggplot2::geom_text(ggplot2::aes(label=CLID), data=centroids, color=thecolors, hjust=0.5,vjust=0.5, size=fs, , fontface="bold")
   }
   g
}


## Loadings plot : Input : data [N,PC] N variables (rows), PC (columns)
## 
## use with associations file: plot of ellipses corresponding to each cluster
## 

#' plotLoadings
#'
#' Plots the two components defined by pc1, pc2 of the matrix of variable loadings coming from a 
#' multivariable analysis, typically a Principal Component Analysis (PCA).
#' It can also plot the ellipses corresponding to each cluster defined by the associations matrix 
#' if not null. (in fact it is the main interest of this function).
#'
#' @param data the matrix of variable loadings coming from a multivariable analysis, typically a Principal Component Analysis (PCA)
#' @param pc1 the fist component of the matrix of variable loadings to be plotted.
#' @param pc2 the second component of the matrix of variable loadings to be plotted.
#' @param EV Eigenvalues vector
#' @param associations the associations matrix that gives for each cluster (column 2) the corresponding buckets (column 1). See \code{getClusters}
#' @param main Change the default plot title on the rigth corner
#' @param onlylabels if TRUE, put only the association names without drawing the cluster contours. Implies that association matrix is provided.
#' @param highlabels if TRUE, put the the association names in blue, and others in grey. Implies that association matrix is provided and fONLYLABELS equal to TRUE.
#' @param gcontour type of contour; possible values are : 'ellipse', 'polygon', 'ellipse2', 'none'
ggplotLoadings <- function (data, pc1=1, pc2=2, EV=NULL, associations=NULL, main="Loadings", onlylabels=FALSE, highlabels=FALSE, gcontour="ellipse" )
{
   P <- data[,c(pc1,pc2)]
   Loadings <- data.frame(IDS=rownames(P), pc1=P[,1], pc2=P[,2])
   xlabs <- colnames(P)[1]; ylabs <- colnames(P)[2]
   if (! is.null(EV)) {
      xlabs <- paste0(xlabs," (",round(EV[pc1],2),"%)"); ylabs <- paste0(ylabs," (",round(EV[pc2],2),"%)")
   }

   fclust <- (!is.null(associations))
   draw.points <- TRUE
   draw.labels <- TRUE
   lw <- 0.3 # linewidth
   ps <- 0.5 # pointsize
   fs <- 4   # fontsize

  if (!fclust || (fclust && onlylabels)) {
      facpc <- 0.7 # threshold value for highlighting loadings
      Loadings$change <- rep("nochange", dim(Loadings)[1])
      Loadings$change[ Loadings$pc2 > facpc*max(Loadings$pc2) ] <- "UP.PC2"
      Loadings$change[ Loadings$pc2 < facpc*min(Loadings$pc2) ] <- "DOWN.PC2"
      Loadings$change[ Loadings$pc1 > facpc*max(Loadings$pc1) ] <- "UP.PC1"
      Loadings$change[ Loadings$pc1 < facpc*min(Loadings$pc1) ] <- "DOWN.PC1"
   }
   IDS <- CLID <- change <- NULL
   g <- ggplot2::ggplot(Loadings, ggplot2::aes(x=pc1, y=pc2)) + ggplot2::ggtitle(main) +
         ggplot2::geom_hline(yintercept=0, color="red", linetype="dashed", , size=0.2) + 
         ggplot2::geom_vline(xintercept=0, color="red", linetype="dashed", , size=0.2) +
         ggplot2::labs(x=xlabs, y=ylabs) + ggplot2::theme_light()

   if (fclust && !onlylabels) {
         g <- g + ggplot2::theme(legend.position="none", text=ggplot2::element_text(size=9), 
              panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
              plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
   }
   if (fclust) {
         if (onlylabels) { # Merged Clusters
              L <- which( rownames(P) %in% associations[,2] )
              Clusters <- data.frame( pc1 = P[L,1], pc2 = P[L,2], CLID = rownames(P)[L] )
              if (highlabels) {
                    g <- g + ggplot2::geom_text(ggplot2::aes(label=IDS),hjust=0.5,vjust=0.5, color="grey", size=3) +
                             ggplot2::geom_text(ggplot2::aes(x=pc1, y=pc2, label=CLID), data=Clusters, hjust=0.5,vjust=0.5, color="blue", size=3)
              } else {
                    g <- g + ggplot2::geom_text(ggplot2::aes(label=IDS, color=change),hjust=0.5,vjust=0.5, size=3) +
                             ggplot2::geom_text(ggplot2::aes(x=pc1, y=pc2, label=CLID, color=change), data=Clusters, hjust=0.5,vjust=0.5, size=3)
              }

         } else { # No Merged Clusters
              L <- which( rownames(P) %in% associations[,1] )
              Clusters <- data.frame( pc1 = P[L,1], pc2 = P[L,2],
                                      CLID = sapply(rownames(P)[L], function(x) { associations[which(associations[,1] == x),2] }) )
              g <- g + ggplot2::geom_text(ggplot2::aes(label=IDS),hjust=0.5,vjust=0.5, color="lightgrey", size=3, check_overlap = TRUE)
              g <- geom_cluster(g,data=Clusters, level=0.8, lw=lw, ps=ps, fs=fs, 
                              draw.contour=gcontour, draw.points=draw.points, draw.labels=draw.labels)
         }
   } else {
         g <- g + ggplot2::geom_text(ggplot2::aes(label=IDS, color=change),hjust=0.5,vjust=0.5, size=3)
   }
   g
}



## Define the 'plot.scores' function for plotting PCA scores.
##  Input : data [N,PC] - N samples (rows), PC (columns)
##          samples - matrix of samples generated by the Rnmr1D function
##          factor - name (string) of the factor ; must be present as a column within the samples matrix
##          level - confidence level for plotting the corresponding ellipse

#' ggplotScores
#'
#' Plots the two components defined by pc1, pc2 of the matrix of scores coming from a 
#' multivariable analysis, typically a Principal Component Analysis (PCA).
#'
#' @param data the matrix of scores coming from a multivariable analysis, typically a Principal Component Analysis (PCA)
#' @param pc1 the fist component of the matrix of variable loadings to be plotted.
#' @param pc2 the second component of the matrix of variable loadings to be plotted.
#' @param groups the vector defining the factorial groups (same dimension as data rows)
#' @param EV Eigenvalues vector
#' @param main the plot main title
#' @param glabels boolean indicating if labels have to be plotted
#' @param psize point size
#' @param gcontour type of contour; possible values are : 'ellipse', 'polygon', 'ellipse2', 'none'
#' @param params parameters depending on the contour type
#' @param colors array of colors
ggplotScores <- function (data, pc1=1, pc2=2, groups=NULL, EV=NULL, main="Scores", glabels=FALSE, psize=3, gcontour="ellipse", params=list(cellipse=0.95), colors=NULL)
{
   S <- data[,c(pc1,pc2)]
   xlabs <- colnames(S)[1]; ylabs <- colnames(S)[2]
   if (! is.null(EV)) {
      xlabs <- paste0(xlabs," (",round(EV[pc1],2),"%)"); ylabs <- paste0(ylabs," (",round(EV[pc2],2),"%)")
   }
   if (is.null(groups)) {
      Scores <- data.frame(IDS=rownames(S), pc1=S[,1], pc2=S[,2])
      g <- ggplot2::ggplot(Scores, ggplot2::aes(pc1, pc2)) + ggplot2::ggtitle(main) +
           ggplot2::geom_point() +
           ggplot2::labs(x=xlabs, y=ylabs, colour = "") + ggplot2::theme_light() + 
           ggplot2::theme(legend.position="none", plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt") )
      IDS <- NULL
      if (glabels)
           g <- g + ggplot2::geom_text(ggplot2::aes(label=IDS),hjust=0.5,vjust=0.5, size=3)
   } else {
      Scores <- data.frame(IDS=rownames(S), pc1=S[,1], pc2=S[,2], groups=groups)
      if (glabels) {
          g <- ggplot2::ggplot(Scores, ggplot2::aes(pc1, pc2, color = groups, fill = groups)) + 
               ggplot2::geom_point() + ggplot2::geom_text(label=Scores$IDS, hjust=0.5,vjust=0.5, size=3)
      } else {
          g <- ggplot2::ggplot(Scores, ggplot2::aes(pc1, pc2, shape=groups, color = groups, fill = groups)) +
               ggplot2::geom_point(size=psize)
          shape_values <- c(15:17, 23,25, 18:22, 0:14, 15:17, 23,25, 18:22, 0:14, 15:17, 23,25, 18:22, 0:14)
          if (nlevels(groups)<=length(shape_values)) {
              g <- g + ggplot2::scale_shape_manual("",values = shape_values )
          }
      }
      if (!is.null(colors)) {
          g <- g + ggplot2::scale_colour_manual(values=colors)
          g <- g + ggplot2::scale_fill_manual(values=colors)
      }
      if (gcontour=="ellipse")  g <- g + ggplot2::stat_ellipse(geom = "polygon", alpha = 0.2, level=params$cellipse, show.legend=FALSE)
      if (gcontour=="ellipse2") g <- g + ggplot2::stat_ellipse(ggplot2::aes(color = groups), level=params$cellipse, show.legend=FALSE)
      if (gcontour=="polygon") {
          .find_hull <- function(Scores) Scores[grDevices::chull(Scores[,2], Scores[,3]), ]
          hulls <- plyr::ddply(Scores, "groups", .find_hull)
          g <- g + ggplot2::geom_polygon(data = hulls, ggplot2::aes(pc1, pc2, fill=factor(groups)), alpha = 0.2, show.legend=FALSE)
      }
      g <- g + ggplot2::labs(x=xlabs, y=ylabs, color = "") + ggplot2::theme_light() + ggplot2::ggtitle(main) + 
               ggplot2::theme(legend.text=ggplot2::element_text(size=9), plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt") )
   }
   g
}

# ggplotEvnorm <- function(evnom, xlab="PC", ylab="% Variance")
# {
#    df<-data.frame(EV=evnom,PC=c(1:length(evnom)))
#    g <- ggplot2::ggplot( df, ggplot2::aes(x=PC, y=EV) ) + ggplot2::geom_bar(stat="identity") + ggplot2::labs(x = xlab, y = ylab) +
#          ggplot2::geom_text(ggplot2::aes(y=0, label=PC), size=4, vjust=1.5) +
#          ggplot2::theme_light() + 
#          ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
#                         panel.grid.minor = ggplot2::element_blank(), 
#                         panel.border = ggplot2::element_rect(color="white"),
#                         axis.ticks.y=ggplot2::element_line(color="black"), 
#                         axis.text.x=ggplot2::element_blank(), 
#                         axis.ticks.x=ggplot2::element_blank())
#    g
# }
# 
# ggplotCompBarplot <- function (M,cp,nbvar,title=NULL)
# {
#     if ((assoc1) && (dim(associations)[1]>1)) {
#        rownames(M)[rownames(M) %in% associations[,1]] <- associations[associations[,1] %in% rownames(M),2]
#     }
#     #r <- round(abs(max(M[,cp])/min(M[,cp])),0)
#     r <- 1
#     pm <- nbvar*r/(1+r)
#     nm <- nbvar - pm
#     N <-M[order(M[,cp],decreasing = T),cp]
#     Np <- length(N[N>=0])
#     Nn <- length(N[N<0])
#     if (Np>=pm) {
#        if (Nn>=nm) { p <- pm;  n <- nm; }        else { n <- Nn; p <- nbvar - n;  }
#     } else {
#        if (Nn>=nm) { p <- Np;  n <- nbvar - p; } else { n <- Nn; p <- Np; }
#     }
#     Vp <- NULL; if (p>0) Vp <- N[1:p] 
#     Vn <- NULL; if (n>0) Vn <- N[(length(N)-n+1):length(N)]
#     V <- c( Vp, Vn )
#     if (is.null(title)) title <- colnames(M)[cp]
#     
#     df <- data.frame( X=names(V), Y=V, group=order(V) )
#     df$group <- factor(df$group, levels = df$group)
#     g <- ggplot2::ggplot(df, ggplot2::aes(x=group, y=Y, fill=group)) + ggplot2::coord_flip() + 
#          ggplot2::labs(x = "", y = "") + ggplot2::ggtitle(colnames(M)[cp]) + 
#          ggplot2::geom_bar(stat = "identity", , show.legend = FALSE) + 
#          ggplot2::geom_text(ggplot2::aes(x=group, y=0, label=X), size=4,  hjust=ifelse(sign(df$Y)>0, 1.1, -0.1)) +
#          ggplot2::theme_light() + theme(axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank(), 
#                                panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
#                                panel.border = ggplot2::element_rect(color="white"), 
#                                axis.line.x=ggplot2::element_line(color="black"), axis.ticks.x=ggplot2::element_line(color="black"))
#     g
# }

#' ggplotPlotly
#'
#' Translate 'ggplot2' graphs to an interactive plotly version
#'
#' @param g The ggplot2 graph object to be translated into an interactive plotly version
#' @param width Width of the plot in pixels (optional, defaults to automatic sizing).
#' @param height Height of the plot in pixels (optional, defaults to automatic sizing)
#' @param textposition Position of the labels on the graphs relative to the points. Possible values are : 'right', 'left', 'top' or 'buttom'
ggplotPlotly <- function(g, width=NULL, height=NULL, textposition = "right")
{
      options(warn=-1)
      gg <- plotly::ggplotly(g, width=width, height=height) %>% plotly::style(textposition = textposition)
      gg$x$config$modeBarButtonsToAdd <- NULL
      gg$x$layout$margin$t <- gg$x$layout$margin$t + 50

      ww <- plotly::as_widget(gg)
      ww$sizingPolicy$padding <- 10
      ww$sizingPolicy$viewer$padding <- 10
      ww$sizingPolicy$browser.padding <- 10
      ww
}
