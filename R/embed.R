# This file is part of EmbedSOM.
#
# Copyright (C) 2018-2020 Mirek Kratochvil <exa.exa@gmail.com>
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

#' Process the cells with SOM into a nice embedding
#' 
#' @param data Data matrix with points that optionally overrides the one from `fsom$data`.
#' @param map Map object in FlowSOM format, to optionally override `fsom$map`
#' @param fsom FlowSOM object with a built SOM (used if data or map are missing)
#' @param smooth Produce smoother (positive values) or more rough approximation (negative values).
#' @param k How many SOM neighbors to take into the whole computation
#' @param adjust How much non-local information to remove (parameter a)
#' @param importance Importance of dimensions that was used to train the SOM
#' @param coords A matrix of embedding-space coordinates that correspond to 'map$codes' (i.e. the "embedded landmarks"). Overrides 'map$grid' if not NULL.
#' @param coordsFn A coordinates-generating function (e.g. 'EmbedSOM::tSNECoords()') that overrides the existing 'map$grid'.
#' @param emcoords Provided for backwards compatibility, will be removed. Either a matrix of embedded coordinates (same number of rows as 'map$codes', and either 2 or 3 columns depending on the SOM grid dimension), or one of 'flat' (default behavior), 'som' (adjust the SOM coords according to U-matrix distances), 'mst' (embed to MST-like structure), 'fsom-mst' (embed to MST that should look exactly like that of FlowSOM), 'tsne' (embed using tSNE from package Rtsne), 'umap' (embed using UMAP from package umap) or 'uwot::umap' (embed using UMAP from package uwot)
#' @param emcoords.pow Provided for backwards compatibility, will be removed. Exaggeration factor (power) of the distances in U-matrix used for some methods of auto-generating emcoords; default 1.
#' @param threads Number of threads used for computation, 0 chooses hardware concurrency, 1 (default) turns off parallelization.
#' @param parallel Boolean flag whether the computation should be parallelized (this flag is just a nice name for 'threads' and does not do anything directly -- default FALSE sets threads=1, TRUE sets threads=0)
#' @return matrix with 2D or 3D coordinates of the embedded cels, depending on the map
#'
#' @examples
#' d <- cbind(rnorm(10000), 3*runif(10000), rexp(10000))
#' colnames(d) <- paste0("col",1:3)
#' map <- EmbedSOM::SOM(d, xdim=10, ydim=10)
#' e <- EmbedSOM::EmbedSOM(data=d, map=map)
#' EmbedSOM::PlotEmbed(e, data=d, 'col1', pch=16)
#' @useDynLib EmbedSOM, .registration = TRUE
#' @export

EmbedSOM <- function(data=NULL, map=NULL, fsom=NULL,
                     smooth=NULL, k=NULL, adjust=NULL,
                     importance=NULL,
                     coordsFn=NULL, coords=NULL,
                     emcoords=NULL, emcoords.pow=1,
                     parallel=F, threads=if (parallel) 0 else 1) {

  if(is.null(map)) {
    if(is.null(fsom)) {
      stop("You need to supply a map.")
    } else map <- fsom$map
  }

  if(is.null(data)) {
    if(is.null(fsom)) {
      stop("You need to supply the data points.")
    } else data <- fsom$data
  }

  if(is.null(coordsFn) && !is.null(emcoords)) {
    if(length(emcoords)!=1)
      coords <- emcoords
    else {
      warn("emcoords parameter will be deprecated, use coordsFn instead.")
      if(emcoords=='flat') coordsFn <- function(x)x
      if(emcoords=='som') coordsFn <- UMatrixCoords(distFn=function(x)x^emcoords.pow)
      if(emcoords=='mst') coordsFn <- MSTCoords(distFn=function(x)x^emcoords.pow)
      if(emcoords=='fsom-mst') coordsFn <- MSTCoords(distFn=function(x)x^emcoords.pow,
        layoutFn=function(coords,...)igraph::layout_with_kk(...))
      if(emcoords=='tsne') coordsFn <- tSNECoords()
      if(emcoords=='umap') coordsFn <- UMAPCoords()
      if(emcoords=='uwot') coordsFn <- uwotCoords()
    }
  }

  if(is.null(coords)) {
    if(!is.null(coordsFn))
      map <- coordsFn(map)
    coords <- map$grid
  }

  if(is.null(coords))
    stop("Missing embedding coordinates (specify coords or coordsFn, or use a map with grid)")

  somdim <- dim(coords)[2]
  if(!(somdim %in% c(2,3)))
    stop("Unsupported embedding dimension (check size of the grid).")

  ncodes <- dim(map$codes)[1]

  ndata <- nrow(data)
  colsUsed <- map$colsUsed
  if(is.null(colsUsed)) colsUsed <- (1:ncol(data))
  dim <- length(colsUsed)

  if(is.null(smooth))
    smooth <- 0

  if(is.null(k)) {
    k <- as.integer(1+sqrt(ncodes))
  }

  if(is.null(adjust)) {
    adjust <- 1
  }

  if (smooth < -3) {
    stop("Value of smooth must be at least -3.")
  }

  boost <- exp(-smooth-1)

  if (k < 2+somdim) {
    stop(paste("Use at least", 2+somdim, "neighbors for sane results!"))
  }

  if (k>ncodes) {
    stop("k must be less than number of codes in SOM!")
  }

  if(adjust<0) {
    stop("adjust must not be negative!");
  }

  if(!is.null(importance)) {
    if(!is.vector(importance) || length(importance)!=length(colsUsed))
      stop("Importance must be null, or a vector that matches colsUsed length")
    points <- t(data[,colsUsed]) * importance
  } else
    points <- t(data[,colsUsed])


  if(somdim==2) {
    if(dim(coords)[1] != ncodes || dim(coords)[2] != 2) {
      stop("Embedding coordinates need to be of dimension (n_codes, 2).")
    }
  } else {
    if(dim(coords)[1] != ncodes || dim(coords)[2] != 3) {
      stop("Embedding coordinates need to be of dimension (n_codes, 3).")
    }
  }

  # validate size of all data matrices that go into C
  if(dim(map$codes)[1] != ncodes || dim(map$codes)[2] != dim)
    stop("wrong size of the codebook")

  # points are already transposed
  if(dim(points)[2] != ndata || dim(points)[1] != dim)
    stop("wrong size of the input data (check out the column names)")

  if(dim(coords)[1] != ncodes || dim(coords)[2] != somdim)
    stop("wrong number of embedding coordinates")

  codes <- t(map$codes)
  coords <- t(coords)

  embedding <- matrix(0, ncol=nrow(data), nrow=somdim)

  res <- .C("C_embedSOM",
    pnthreads=as.integer(threads),

    psomdim=as.integer(somdim),
    pn=as.integer(ndata),
    pncodes=as.integer(ncodes),
    pdim=as.integer(dim),

    pdistf=as.integer(map$distf),

    pboost=as.single(boost),
    pneighbors=as.integer(k),
    padjust=as.single(adjust),

    points=as.single(points),
    koho=as.single(codes),
    coords=as.single(coords),

    embedding=as.single(embedding))

  matrix(res$embedding,
    byrow=T,
    nrow=nrow(data),
    ncol=somdim,
    dimnames=list(NULL,paste0('EmbedSOM',seq_len(somdim))))
}
