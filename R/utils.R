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
#' @export
NormalizeColor <- function(data, low=0, high=1, pow=0) {
  ps <- quantile(data[!is.nan(data)], c(low,high))
  ps <- pmin(pmax((data-ps[1])/(ps[2]-ps[1]), 0),1)^(2^pow)
  ps[is.nan(ps)] <- 0
  ps
}

#' Marker expression palette generator defined from ColorBrewer's RdYlBu
#' @export
ExpressionPalette <- function(n, alpha=1) {
  # proudly taken from colorBrewer (thanks!)
  pal <- rev(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))

  grDevices::adjustcolor(alpha=alpha,
    col=grDevices::colorRampPalette(pal, space='Lab', interpolate='linear')(n))
}

#' An acceptable cluster color palette
#' @export
ClusterPalette <- function(n, vcycle=c(1,0.7), scycle=c(0.7,1), alpha=1)
{
  hsv(alpha=alpha, h=c(0:(n-1))/n, v=vcycle, s=scycle)
}

#' Identity on whatever
PlotId <- function(x){x}

#' Helper function for plotting the embedding
#' 
#' Takes the 'embed' object which is the output of EmbedSOM, together with a
#' multitude of arguments that override how the plotting is done.
#'
#' @param embed The embedding from EmbedSOM
#' @param value The column of data to use for plotting the value
#' @param red,green,blue The same for RGB components
#' @param fv,fr,fg,fb Functions to transform the values before they are normalized
#' @param powv,powr,powg,powb Adjustments of the value plotting
#' @param limit Low/high offset for NormalizeColor
#' @param clust,nclust integer cluster label column/data, optional cluster number
#' @param alpha Default alpha value
#' @param col Different coloring, if supplied
#' @param pch,cex Parameters for point plots
#' @param exlim Extra border around the embedding, default 1
#' @param data Data matrix, taken from fsom parameter by default
#' @param fsom FlowSOM object
#' @param xdim,ydim Sizes of the SOM, taken from fsom by default
#' @export
PlotEmbed <- function(embed,
  value=0, red=0, green=0, blue=0,
  fr=PlotId, fg=PlotId, fb=PlotId, fv=PlotId,
  powr=0, powg=0, powb=0, powv=0,
  clust=NULL, nclust=0,
  nbin=256, maxDens=NULL,
  limit=0.01, pch='.', alpha=NULL, cex=1, exlim=1, fsom, data, xdim, ydim, col) {
  if(missing(data)) {
    data <- fsom$data
  }
  if(missing(xdim)) {
    xdim <- fsom$map$xdim
  }
  if(missing(ydim)) {
    ydim <- fsom$map$ydim
  }
  if(missing(col)) {
    if (!is.null(clust)) {
      if(length(clust)==1) cdata <- data[,clust]
      else cdata <- clust
      if(nclust==0) nclust <- {
        tmp <- cdata
        tmp[is.nan(tmp)]<-0
        max(tmp)
      }
      cdata[cdata<1 | cdata>nclust] <- NaN #produce NaNs instead of skipping
      col <- ClusterPalette(nclust, alpha=alpha)[cdata]
    } else if(value==0 & red==0 & green==0 & blue==0) {
      if(is.null(alpha)) alpha <- 1
      xbin <- cut(embed[,1], -exlim+(xdim+2*exlim)*c(0:nbin)/nbin, labels=F)
      ybin <- cut(embed[,2], -exlim+(ydim+2*exlim)*c(0:nbin)/nbin, labels=F)

      dens <- tabulate(xbin+(nbin+1)*ybin)[xbin+(nbin+1)*ybin]
      if(!is.null(maxDens)) dens[dens>maxDens] <- maxDens
      dens <- log(dens+1)
      pal <- cut(dens, length(dens), labels=F)
      n <- length(dens)
      col <- ExpressionPalette(256, alpha=alpha)[1+as.integer(255*pal/n)]
    } else if(value==0) {
      if(is.null(alpha)) alpha <- 0.5
      col <- rgb(
        if(red>0)   NormalizeColor(fr(data[,red]),   limit, 1-limit, powr) else 0,
        if(green>0) NormalizeColor(fg(data[,green]), limit, 1-limit, powg) else 0,
        if(blue>0)  NormalizeColor(fb(data[,blue]),  limit, 1-limit, powb) else 0,
      alpha)
    } else {
      if(is.null(alpha)) alpha <- 0.5
      col <- ExpressionPalette(256,alpha=alpha)[1+255*NormalizeColor(fv(data[,value]), limit, 1-limit, powv)]
    }
  }

  plot(
    embed,
    cex=cex,
    pch=pch,
    col=col,
    xaxt='n',
    yaxt='n',
    xaxs='i',
    yaxs='i',
    xlim=c(-exlim, xdim+exlim-1),
    ylim=c(-exlim, ydim+exlim-1));
}
