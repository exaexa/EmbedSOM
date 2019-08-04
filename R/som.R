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
#
# This particular file is based on FlowSOM (C) Sofie Van Gassen (2016-)


#' Build a self-organizing map
#'
#' @param data  Matrix containing the training data
#' @param xdim  Width of the grid
#' @param ydim  Hight of the grid
#' @param zdim  Depth of the grid, causes grid to be 3D
#' @param rlen  Number of times to loop over the training data for each MST
#' @param alphaA Start and end learning rate
#' @param alphaB Start and end learning rate for the second radius
#' @param radiusA Start and end radius
#' @param radiusB Start and end radius (make sure it's larger than radiusA)
#' @param negRadius easy way to set radiusB as a multiple of default radius
#'                  (use lower value for higher dimensions)
#' @param negAlpha the same for alphaB
#' @param init  Initialize cluster centers in a non-random way
#' @param initf Use the given initialization function if init==T
#'              (default: Initialize_PCA)
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev,
#'              4=cosine)
#' @param codes Cluster centers to start with
#' @param importance array with numeric values. Parameters will be scaled
#'                   according to importance
#' @param nhbr.method Way of computing grid distances, passed as method= to dist() function. Default 'maximum' (square neighborhoods); use 'euclidean' for round neighborhoods.
#' @param noMapping If true, do not produce mapping (default F). Useful for online/streaming use.
#'
#' @return A map, which is a list containing all parameter settings and results
#'
#' @seealso FlowSOM::SOM
#'
#' @references This code is strongly based on the \code{kohonen} package.
#'             R. Wehrens and L.M.C. Buydens, Self- and Super-organising Maps
#'             in R: the kohonen package J. Stat. Softw., 21(5), 2007
#' @useDynLib EmbedSOM, .registration = TRUE
#' @export

SOM <- function (data, xdim=10, ydim=10, zdim=NULL, rlen=10,
    alphaA=c(0.05, 0.01), radiusA = stats::quantile(nhbrdist, 0.67) * c(1, 0),
    alphaB=alphaA*c(-negAlpha,-0.01*negAlpha), radiusB = negRadius*radiusA,
    init=FALSE, initf=Initialize_PCA, distf=2,
    codes=NULL, importance = NULL, nhbr.method='maximum',
    negRadius=1.33, negAlpha=0.1,
    noMapping=F) {

    somdim <- 2
    if(!is.null(zdim)) somdim <- 3

    if (!is.null(codes)){
      if((ncol(codes) != ncol(data)) | !(
        (somdim==2 && (nrow(codes) == xdim * ydim)) |
        (somdim==3 && (nrow(codes) == xdim * ydim * zdim))
      )){
        stop("If codes is not NULL, it should have the same number of columns
             as the data and the number of rows should correspond with
             xdim*ydim (*zdim)")
      }
    }

    if(!is.null(importance)){
        data <- data * rep(importance,each=nrow(data))
    }

    # Initialize the grid
    if(somdim==2) grid <- expand.grid(seq_len(xdim),seq_len(ydim))
    else if(somdim==3) grid <- expand.grid(seq_len(xdim), seq_len(ydim), seq_len(zdim))
    else stop("somdim must be either 2 or 3")

    nCodes <- nrow(grid)

    if(is.null(codes)){
        if(init){
            codes <- initf(data, xdim, ydim, zdim)
            message("Initialization ready\n")
        } else {
            codes <- data[sample(1:nrow(data), nCodes, replace = FALSE), , drop = FALSE]
        }
    }

    # Initialize the neighbourhood
    nhbrdist <- as.matrix(stats::dist(grid, method = nhbr.method))

    # Compute the SOM
    res <- .C("es_C_SOM",
        data = as.single(data),
        codes = as.single(codes),
        nhbrdist = as.single(nhbrdist),
        alphaA = as.single(alphaA),
        radiusA = as.single(radiusA),
        alphaB = as.single(alphaB),
        radiusB = as.single(radiusB),
        n = as.integer(nrow(data)),
        px = as.integer(ncol(data)),
        ncodes = as.integer(nCodes),
        rlen = as.integer(rlen),
        distf = as.integer(distf))

    codes <- matrix(res$codes, nrow(codes), ncol(codes))
    colnames(codes) <- colnames(data)

    if(noMapping) mapping <- NULL
    else mapping <- MapDataToCodes(codes,data)

    list(xdim=xdim, ydim=ydim, zdim=zdim, rlen=rlen,
        alphaA=alphaA, radiusA=radiusA,
        alphaB=alphaB, radiusB=radiusB,
        init=init, distf=distf,
        nhbr.method=nhbr.method,
        grid=grid, codes=codes, mapping=mapping, nNodes=nCodes)
}



#' Assign nearest node to each datapoint
#
#' @param codes matrix with nodes of the SOM
#' @param newdata datapoints to assign
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev,
#'              4=cosine)
#'
#' @return Array with nearest node id for each datapoint
#' @export
MapDataToCodes <- function (codes, newdata, distf=2) {

    nnCodes <- .C("es_C_mapDataToCodes",
        as.single(newdata[,colnames(codes)]),
        as.single(codes),
        as.integer(nrow(codes)),
        as.integer(nrow(newdata)),
        as.integer(ncol(codes)),
        nnCodes = integer(nrow(newdata)),
        nnDists = double(nrow(newdata)),
        distf = as.integer(distf))
    cbind(nnCodes$nnCodes, nnCodes$nnDists)
}

#' Create a grid from first 2 PCA components
#'
#' @param   data matrix in which each row represents a point
#' @param   xdim,ydim,zdim Dimensions of the SOM grid
#' @return  array containing the selected selected rows
#' @export
Initialize_PCA <- function(data, xdim, ydim, zdim=NULL){
    if(!is.null(zdim)) stop("No support for 3D PCA yet")

    pca <- stats::prcomp(data, rank.=2, retx=F)
    sdev_scale <- 5 # scale out to 5-times standard deviation, which should cover the data nicely
    ax1 <- t(matrix(pca$rotation[,1] * sdev_scale * pca$sdev,
             nrow=ncol(data),
             ncol=xdim * ydim)) *
             (2 * rep(c(1:xdim) - 1, times=ydim) / (xdim - 1) - 1)
    ax2 <- t(matrix(pca$rotation[,2] * sdev_scale * pca$sdev,
             nrow=ncol(data),
             ncol=xdim * ydim)) *
             (2 * rep(c(1:ydim) - 1, each=xdim) / (ydim - 1) - 1)

    t(matrix(pca$center, nrow=ncol(data), ncol=xdim * ydim)) + ax1 + ax2
}

