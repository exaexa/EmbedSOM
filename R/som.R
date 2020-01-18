# This file is part of EmbedSOM.
#
# Copyright (C) 2018-2020 Mirek Kratochvil <exa.exa@gmail.com>
#
# A part of the functionality is based on FlowSOM,
# Copyright (C) 2016-2019 Sofie Van Gassen et al.
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


#' Build a self-organizing map
#'
#' @param data  Matrix containing the training data
#' @param xdim  Width of the grid
#' @param ydim  Hight of the grid
#' @param zdim  Depth of the grid, causes grid to be 3D
#' @param batch Use batch training (default FALSE chooses online training, more compatible with FlowSOM)
#' @param rlen  Number of training epochs; or number of times to loop over the training data in online training
#' @param alphaA Start and end learning rate for online learning (only for online training)
#' @param alphaB Start and end learning rate for the second radius (only for online training)
#' @param radiusA Start and end radius
#' @param radiusB Start and end radius (only for online training; make sure it's larger than radiusA)
#' @param negRadius easy way to set radiusB as a multiple of default radius
#'                  (use lower value for higher dimensions)
#' @param negAlpha the same for alphaB
#' @param epochRadii Vector of length 'rlen' with precise epoch radii (only for batch training)
#' @param init  Initialize cluster centers in a non-random way
#' @param initf Use the given initialization function if init==T
#'              (default: Initialize_PCA)
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev)
#' @param codes Cluster centers to start with
#' @param importance array with numeric values. Parameters will be scaled
#'                   according to importance
#' @param coordsFn Function to generate/transform grid coordinates (e.g. 'EmbedSOM::tSNECoords()'). If NULL (default), the grid is the canonical SOM grid.
#' @param nhbr.method Way of computing grid distances, passed as method= to dist() function. Default 'maximum' (square neighborhoods); use 'euclidean' for round neighborhoods.
#' @param noMapping If true, do not produce mapping (default F). Useful for online/streaming use.
#' @param threads Number of threads of batch training (has no effect on online training). Defaults to 0 (chooses maximum available hardware threads) if parallel=TRUE or 1 (single thread) if parallel=FALSE. Passed to 'MapDataToCodes'.
#' @param parallel Parallelize batch training (has no effect on online training) by setting appropriate 'threads'. Defaults to FALSE. Use 'batch=T' for fully parallelized version, online training is not parallelizable. Passed to 'MapDataToCodes'.
#' @return A map useful for embedding ('EmbedSOM' function) or further analysis.
#'
#' @seealso FlowSOM::SOM
#'
#' @useDynLib EmbedSOM, .registration = TRUE
#' @export

SOM <- function (data, xdim=10, ydim=10, zdim=NULL, batch=F, rlen=10,
    alphaA=c(0.05, 0.01), radiusA=stats::quantile(nhbrdist, 0.67) * c(1, 0),
    alphaB=alphaA*c(-negAlpha,-0.1*negAlpha), radiusB = negRadius*radiusA,
    negRadius=1.33, negAlpha=0.1,
    epochRadii=seq(radiusA[1], radiusA[2], length.out=rlen),
    init=FALSE, initf=Initialize_PCA, distf=2,
    codes=NULL, importance = NULL, coordsFn = NULL,
    nhbr.method='maximum',
    noMapping=F, parallel=F, threads=if (parallel) 0 else 1) {

    if(threads!=1 && batch==F)
      warning("The used SOM training version is not parallelizable. Perhaps you want to use batch=T?")

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

    if(is.null(colnames(data)))
      stop("Incoming data must have correctly named columns")

    if(!is.null(importance)){
      if(!is.vector(importance) || length(importance)!=ncol(data))
        stop("Importance must be null, or a vector of the column-size of data.")
      points <- t(data) * importance
    } else points <- t(data)

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

    # validate size of data that go into C
    if(dim(codes)[1] != nCodes || dim(codes)[2] != ncol(data))
      stop("wrong size of SOM codebook (check column names of data)")

    if(!is.numeric(epochRadii) || length(epochRadii)<rlen)
      stop("epochRadii must be a numeric vector of length rlen")

    codes <- t(codes)

    # Compute the SOM
    if(batch) res <- .C("es_C_BatchSOM",
      nthreads=as.integer(threads),
      data = as.single(points),
      codes = as.single(codes),
      nhbrdist = as.single(nhbrdist),
      epochRadii = as.single(epochRadii),
      n = as.integer(nrow(data)),
      dim = as.integer(ncol(data)),
      kohos = as.integer(nCodes),
      rlen = as.integer(rlen),
      distf = as.integer(distf))
    else res <- .C("es_C_SOM",
      data = as.single(points),
      codes = as.single(codes),
      nhbrdist = as.single(nhbrdist),
      alphaA = as.single(alphaA),
      radiusA = as.single(radiusA),
      alphaB = as.single(alphaB),
      radiusB = as.single(radiusB),
      n = as.integer(nrow(data)),
      dim = as.integer(ncol(data)),
      kohos = as.integer(nCodes),
      rlen = as.integer(rlen),
      distf = as.integer(distf))

    codes <- matrix(res$codes, byrow=T, nrow=ncol(codes), ncol=nrow(codes))
    colnames(codes) <- colnames(data)

    if(noMapping) mapping <- NULL
    else mapping <- MapDataToCodes(codes,data, parallel=parallel, threads=threads)

    (if(!is.null(coordsFn)) coordsFn else function(x)x)(list(
      xdim=xdim, ydim=ydim, zdim=zdim, rlen=rlen,
      alphaA=alphaA, radiusA=radiusA,
      alphaB=alphaB, radiusB=radiusB,
      init=init, distf=distf,
      nhbr.method=nhbr.method,
      grid=grid, codes=codes, mapping=mapping, nNodes=nCodes))
}

#' Train a Growing Quadtree Self-Organizing Map
#'
#' @param data Input data matrix
#' @param init.dim Initial size of the SOM, default c(3,3)
#' @param target_codes Make the SOM grow linearly to at most this amount of nodes (default 100)
#' @param rlen Number of training iterations
#' @param radius Start and end training radius, as in 'SOM'
#' @param epochRadii Precise radii for each epoch (must be of length 'rlen')
#' @param coords Tree coordinates of the initial SOM nodes.
#' @param codes Initial codebook
#' @param coordsFn Function to generate/transform grid coordinates (e.g. 'EmbedSOM::tSNECoords()'). If NULL (default), the grid is the grid generated by GQTSOM.
#' @param importance Weights of input data dimensions
#' @param distf Distance measure to use in input data space (1=manhattan, 2=euclidean, 3=chebyshev)
#' @param nhbr.distf Distance measure to use in output space (as in 'distf')
#' @param noMapping If TRUE, do not compute the assignment of input data to SOM nodes
#' @param threads Number of threads to use for training. Defaults to 0 (chooses maximum available hardware threads) if parallel=TRUE or 1 (single thread) if parallel=FALSE.
#' @param parallel Parallelize the training by setting appropriate 'threads'. Defaults to FALSE.
#' @export
GQTSOM <- function(data, init.dim=c(3,3), target_codes=100, rlen=10,
  radius=c(sqrt(sum(init.dim^2)),0.5), epochRadii=seq(radius[1], radius[2], length.out=rlen),
  coords=NULL, codes=NULL, coordsFn = NULL, importance=NULL,
  distf=2, nhbr.distf=2,
  noMapping=F, parallel=F, threads=if (parallel) 0 else 1) {

  xdim <- NULL
  ydim <- NULL

  if(is.null(coords)) {
    if(length(init.dim)!=2)
      stop("Either coords must be initialized, or init.dim must be a vector of 2 integers")
    coords <- as.matrix(cbind(level=as.integer(0),
                              expand.grid(x=1:init.dim[1],
                                          y=1:init.dim[2])))
    xdim <- init.dim[1]
    ydim <- init.dim[2]
  } else if(dim(coords)[1] < 2 || dim(coords)[2]!=3) {
    stop("invalid specified coords! Specify at least 2 coords in a 3-column matrix.")
  } else {
    coords <- as.integer(coords)
    xdim <- dim(coords)[1]
    ydim <- 1
  }

  if(any(coords[,1] < 0))
    stop("Coordinate levels must not be negative!")

  if(is.null(codes)) {
    ncodes <- nrow(coords)
    codes <- data[sample(dim(data)[1], ncodes),,drop=F]
  }

  if(dim(codes)[2]!=dim(data)[2])
    stop("Codes and data must have the same dimension!")

  if(dim(codes)[1]!=dim(coords)[1])
    stop("Different number of codes and coords specified!")

  points <- t(data) * (if(is.null(importance)) 1
                       else if(length(importance)!=dim(codes)[2])
                         stop("Importance must be either null or a vector of the column-size of data")
                         else importance)

  if(target_codes<dim(codes)[1])
    stop("target_codes too low!")

  out.kohos <- as.integer(target_codes)
  out.codes <- single(target_codes*dim(codes)[2])
  out.coords <- integer(target_codes*3)
  out.emcoords <- single(target_codes*2)
  codes <- t(codes)
  coords <- t(coords)

  res <- .C("es_C_GQTSOM",
    nthreads=as.integer(threads),
    data=as.single(points),
    coords=as.integer(coords),
    codes=as.single(codes),
    radii=as.single(epochRadii),
    out.kohos=as.integer(out.kohos),
    out.codes=as.single(out.codes),
    out.coords=as.integer(out.coords),
    out.emcoords=as.single(out.emcoords),
    n=as.integer(ncol(points)),
    dim=as.integer(nrow(codes)),
    kohos=as.integer(ncol(codes)),
    rlen=as.integer(rlen),
    distf=as.integer(distf),
    nhbr.distf=as.integer(nhbr.distf))

  codes <- matrix(res$out.codes[1:(res$out.kohos*nrow(codes))],
                  byrow=T, nrow=res$out.kohos, ncol=nrow(codes))
  colnames(codes) <- colnames(data)

  if(noMapping) mapping <- NULL
  else mapping <- MapDataToCodes(codes, data, parallel=parallel, threads=threads)

  (if(!is.null(coordsFn)) coordsFn else function(x)x)(list(
    xdim=xdim, ydim=ydim, rlen=rlen,
    distf=distf, nhbr.distf=nhbr.distf,
    epochRadii=epochRadii,
    grid=matrix(res$out.emcoords[1:(res$out.kohos*2)],
                byrow=T, nrow=res$out.kohos, ncol=2),
    coords=matrix(res$out.coords[1:(res$out.kohos*3)],
                  byrow=T, nrow=res$out.kohos, ncol=3),
    codes=codes,
    mapping=mapping,
    nNodes=res$out.kohos))
}

#' Assign nearest node to each datapoint
#'
#' @param codes matrix with nodes of the SOM
#' @param data datapoints to assign
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev,
#'              4=cosine)
#' @param threads Number of threads used for computation, 0 chooses hardware concurrency, 1 (default) turns off parallelization.
#' @param parallel Boolean flag whether the computation should be parallelized (this flag is just a nice name for 'threads' and does not do anything directly -- default FALSE sets threads=1, TRUE sets threads=0)
#' @return array with nearest node id for each datapoint
#'
#' @export
MapDataToCodes <- function (codes, data, distf=2, parallel=F, threads=if (parallel) 0 else 1) {

    colsToUse <- colnames(codes)
    if(is.null(colsToUse))
      stop("invalid codebook column names")
    if(!all(colsToUse %in% colnames(data)))
      stop("codebook does not match data")

    data <- t(data[,colsToUse])
    codes <- t(codes)

    nnCodes <- .C("es_C_mapDataToCodes",
        threads = as.integer(threads),
        points = as.single(data),
        koho = as.single(codes),
        pn = as.integer(ncol(data)),
        pdim = as.integer(nrow(data)),
        pkohos = as.integer(ncol(codes)),
        mapping = integer(ncol(data)),
        dists = single(ncol(data)),
        distf = as.integer(distf))
    cbind(1+nnCodes$mapping, nnCodes$dists)
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

