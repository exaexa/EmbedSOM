# This file is part of EmbedSOM.
#
# EmbedSOM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EmbedSOM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with EmbedSOM. If not, see <https://www.gnu.org/licenses/>.


#' Helper for computing colors for embedding plots
#'
#' @param data Vector of scalar values to normalize between 0 and 1
#' @param low,high Quantiles to be clamped to 0 and 1 (try low=0.05, high=0.95)
#' @param pow The scaled data are transformed to data^(2^pow). If set to 0,
#'            nothing happens. Positive values highlight differences in the
#'            data closer to 1, negative values highlight differences closer to 0.
#' @export
NormalizeColor <- function(data, low=0, high=1, pow=0) {
  ps <- stats::quantile(data[!is.nan(data)], c(low,high))
  ps <- pmin(pmax((data-ps[1])/(ps[2]-ps[1]), 0),1)^(2^pow)
  ps[is.nan(ps)] <- 0
  ps
}

#' Marker expression palette generator based off ColorBrewer's RdYlBu,
#' only better for plotting of half-transparent cells
#'
#' @param n How many colors to generate
#' @param alpha Opacity of the colors
#' @export
ExpressionPalette <- function(n, alpha=1) {
  pal <- rev(c(
    "#A50026",
    "#D73027",
    "#F46D43",
    "#FDAE61",
    "#FEE090",
    "#FFFFA8", # this was darkened a bit
    "#B8D9C8", # this was converted to gray from almost-white
    "#91C3E2", # and the rest got darkened a bit
    "#649DD1",
    "#3565B4",
    "#212695"))

  grDevices::adjustcolor(alpha=alpha,
    col=grDevices::colorRampPalette(pal, space='Lab', interpolate='linear')(n))
}

#' An acceptable cluster color palette
#'
#' @param n How many colors to generate
#' @param vcycle,scycle Small vectors with cycles of saturation/value for hsv
#' @param alpha Opacity of the colors
#' @export
ClusterPalette <- function(n, vcycle=c(1,0.7), scycle=c(0.7,1), alpha=1)
{
  grDevices::hsv(alpha=alpha, h=c(0:(n-1))/n, v=vcycle, s=scycle)
}

#' Identity on whatever
#'
#' @param x Just the x.
#' @return The x.
PlotId <- function(x){x}

#' Helper function for plotting the embedding
#'
#' Takes the 'embed' object which is the output of EmbedSOM, together with a
#' multitude of arguments that override how the plotting is done.
#'
#' @param embed The embedding from EmbedSOM
#' @param data Data matrix, taken from fsom parameter by default
#' @param fsom FlowSOM object
#' @param value The column of data to use for plotting the value
#' @param red,green,blue The same for RGB components
#' @param fv,fr,fg,fb Functions to transform the values before they are normalized
#' @param powv,powr,powg,powb Adjustments of the value plotting
#' @param nbin,maxDens,fdens Parameters of density calculation, see PlotData
#' @param limit Low/high offset for NormalizeColor
#' @param clust,nclust integer cluster label column/data, optional cluster number
#' @param alpha Default alpha value
#' @param col Different coloring, if supplied
#' @param pch,cex Parameters for point plots
#' @param cluster.colors Function to generate cluster colors, default ClusterPalette
#' @param expression.colors Function to generate expression color scale, default ExpressionPalette
#' @param ... Extra params passed to plot(...)
#' @export
PlotEmbed <- function(embed,
  value=0, red=0, green=0, blue=0,
  fr=PlotId, fg=PlotId, fb=PlotId, fv=PlotId,
  powr=0, powg=0, powb=0, powv=0,
  clust=NULL, nclust=0,
  nbin=256, maxDens=NULL, fdens=sqrt,
  limit=0.01, pch='.', alpha=NULL, cex=1, fsom, data, col,
  cluster.colors=ClusterPalette,
  expression.colors=ExpressionPalette, ...) {
  if(missing(col)) {
    if(dim(embed)[2]!=2) stop ("PlotEmbed only works for 2-dimensional embedding")

    if (!is.null(clust)) {
      if(length(clust)==1) {
        if(missing(data)) {
          data <- fsom$data
        }
        cdata <- data[,clust]
      }
      else cdata <- clust
      if(nclust==0) nclust <- {
        tmp <- cdata
        tmp[is.nan(tmp)]<-0
        max(tmp)
      }
      cdata[cdata<1 | cdata>nclust] <- NaN #produce NaNs instead of skipping
      col <- cluster.colors(nclust, alpha=alpha)[cdata]
    } else if(value==0 & red==0 & green==0 & blue==0) {
      if(is.null(alpha)) alpha <- 1
      mins <- apply(embed,2,min)
      maxs <- apply(embed,2,max)
      xbin <- cut(embed[,1], mins[1]+(maxs[1]-mins[1])*c(0:nbin)/nbin, labels=F)
      ybin <- cut(embed[,2], mins[2]+(maxs[2]-mins[2])*c(0:nbin)/nbin, labels=F)

      dens <- tabulate(xbin+(nbin+1)*ybin)[xbin+(nbin+1)*ybin]
      if(!is.null(maxDens)) dens[dens>maxDens] <- maxDens
      dens <- fdens(dens)
      pal <- cut(dens, length(dens), labels=F)
      n <- length(dens)
      col <- expression.colors(256, alpha=alpha)[1+as.integer(255*pal/n)]
    } else if(value==0) {
      if(missing(data)) {
        data <- fsom$data
      }
      if(is.null(alpha)) alpha <- 0.5
      col <- grDevices::rgb(
        if(red>0)   NormalizeColor(fr(data[,red]),   limit, 1-limit, powr) else 0,
        if(green>0) NormalizeColor(fg(data[,green]), limit, 1-limit, powg) else 0,
        if(blue>0)  NormalizeColor(fb(data[,blue]),  limit, 1-limit, powb) else 0,
      alpha)
    } else {
      if(missing(data)) {
        data <- fsom$data
      }
      if(is.null(alpha)) alpha <- 0.5
      col <- expression.colors(256,alpha=alpha)[1+255*NormalizeColor(fv(data[,value]), limit, 1-limit, powv)]
    }
  }

  graphics::plot(
    embed,
    cex=cex,
    pch=pch,
    col=col,
    xaxt='n',
    yaxt='n',
    ...)
}

#' Export a data frame for plotting with marker intensities and density.
#'
#' @param embed,fsom,data,cols The embedding data, columns to select
#' @param names Column names for output
#' @param normalize List of columns to normalize using NormalizeColor, default all
#' @param qlimit,qlow,qhigh,pow Parameters for the normalization
#' @param vf Custom value-transforming function
#' @param density Name of the density column
#' @param densBins Number of bins for density calculation
#' @param densLimit Upper limit of density (prevents outliers)
#' @param fdens Density-transforming function; default sqrt
#' @export
PlotData <- function(embed,
  fsom, data=fsom$data, cols, names,
  normalize=cols, qlimit=0, qlow=qlimit, qhigh=1-qlimit, pow=0, vf=PlotId,
  density='Density', densBins=256, densLimit=NULL, fdens=sqrt
  ) {
  if(dim(embed)[2]!=2) stop ("PlotData only works for 2-dimensional embedding")

  if(missing(cols)) {
    cols <- colnames(data)
  }

  df <- data.frame(EmbedSOM1=embed[,1], EmbedSOM2=embed[,2])

  if(is.null(cols) || cols==FALSE) {
    #no cols to add :]
  } else {
    ddf <- data.frame(data[,cols])
    if(missing(names)) {
      if(missing(fsom)) names <- colnames(data)[cols]
      else names <- fsom$prettyColnames[cols]
    }

    colnames(ddf) <- cols
    cols <- colnames(ddf) # you may feel offended but I'm ok. :-/

    ncol <- length(normalize)
    qlow <- rep_len(qlow, ncol)
    qhigh <- rep_len(qhigh, ncol)
    pow <- rep_len(pow, ncol)
    vf <- rep_len(c(vf), ncol)

    for(i in c(1:length(normalize)))
      ddf[,normalize[i]] <- NormalizeColor(
        vf[[i]](ddf[,normalize[i]]),
        qlow[i], qhigh[i], pow[i])

    colnames(ddf) <- names
    df <- data.frame(df, ddf)
  }

  if(!is.null(density)) {
    densBins <- rep_len(densBins, 2)
    xbin <- cut(embed[,1], breaks=densBins[1], labels=F)
    ybin <- cut(embed[,2], breaks=densBins[2], labels=F)

    dens <- tabulate(xbin+(densBins[1]+1)*ybin)[xbin+(densBins[1]+1)*ybin]
    if(!is.null(densLimit)) dens[dens>densLimit] <- densLimit
    n <- length(dens)
    densf <- data.frame(density=cut(fdens(dens), n, labels=F))
    colnames(densf)[1]<-density
    df <- data.frame(df, densf)
  }

  df
}

# this is required because of ggplot aesthetics syntax isn't properly
# recognized by `R CMD check`.
utils::suppressForeignCheck(c("EmbedSOM1", "EmbedSOM2"))

#' Wrap PlotData result in ggplot object.
#'
#' This creates a ggplot2 object for plotting.
#' Use:
#'
#' PlotGG(...) + geom_point()
#'
#' Slight point style modification is recommended:
#'
#' PlotGG(...) + geom_point(aes(color=yourColName), alpha=.3, size=.3)
#'
#' @param embed,fsom Embedding data
#' @param ... Extra arguments passed to PlotData
#' @export
PlotGG <- function(embed, fsom, ...) {
  ggplot2::ggplot(PlotData(embed, fsom, ...)) +
    ggplot2::aes(EmbedSOM1, EmbedSOM2)
}

#' The ggplot2 scale gradient from ExpressionPalette.
#'
#' @example EmbedSOM::PlotGG(...) + EmbedSOM::ExpressionGradient(guide=F)
#'
#' @param ... Arguments passed to ggplot2::scale_color_gradientn
#' @export
ExpressionGradient <- function(...) {
	ggplot2::scale_color_gradientn(colors=ExpressionPalette(256), ...)
}
