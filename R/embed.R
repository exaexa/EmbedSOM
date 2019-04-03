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

#' Process the cells with SOM into a nice embedding
#' 
#' @param fsom FlowSom object with a built SOM
#' @param data Data matrix with points that optionally overrides the one from `fsom$data`
#' @param map Map object in FlowSOM format, to optionally override `fsom$map`
#' @param smooth Produce smoother (positive values) or more rough approximation (negative values).
#' @param k How many SOM neighbors to take into the whole computation
#' @param adjust How much non-local information to remove (parameter a)
#' @param importance Importance of dimensions that was used to train the SOM
#' @param emcoords Either a matrix of embedded coordinates (same number of rows as map$coords, and either 2 or 3 columns depending on the SOM grid dimension), or one of 'flat' (default behavior), 'som' (adjust the SOM coords according to U-matrix distances), 'mst' (embed to MST-like structure), 'fsom-mst' (embed to MST that should look exactly like that of FlowSOM), 'tsne' (embed using tSNE from package Rtsne), 'umap' (embed using UMAP from package umap) or 'uwot::umap' (embed using UMAP from package uwot)
#' @param emcoords.pow Exaggeration factor (power) of the distances in U-matrix used for some methods of auto-generating emcoords; default 1.
#' @return matrix with 2D or 3D coordinates of the embedded cels, depending on the map
#'
#' @useDynLib EmbedSOM, .registration = TRUE
#' @export

EmbedSOM <- function(fsom=NULL, smooth=NULL, k=NULL, adjust=NULL,
                     data=NULL, map=NULL, importance=NULL,
                     emcoords='flat', emcoords.pow=1) {
  #TODO: validate the sizes of data, colsUsed and codes.

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

  somdim <- 2
  if(!is.null(map$zdim)) somdim <- 3

  x <- map$xdim
  y <- map$ydim
  z <- map$zdim

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

  if (smooth< -30) {
    stop("Value of smooth must be at least -30.")
  }

  boost <- ((1+sqrt(5))/2)^(smooth-2)

  if (k<1+somdim) {
    stop(paste("Use at least",somdim+1,"neighbors for sane results!"))
  }

  if (k>ncodes) {
    stop("k must be less than number of codes in SOM!")
  }

  if(adjust<0) {
    stop("adjust must not be negative!");
  }

  if(!is.null(importance))
    points <- t(data[,colsUsed] * rep(importance, each=nrow(data)))
  else
    points <- t(data[,colsUsed])

  if(length(emcoords)==1) {
    theGrid <- as.matrix(map$grid)

    if(emcoords=='flat') {
      emcoords <- theGrid
    } else if(emcoords=='som') {
      emcoords <- igraph::layout_with_kk(coords=theGrid, dim=somdim,
        igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
          as.matrix(igraph::as_adjacency_matrix(
            igraph::make_lattice(if(somdim==2) {c(x,y)} else {c(x,y,z)})))
          *
          as.matrix(dist(map$codes))^emcoords.pow))
    } else if(emcoords=='mst') {
      emcoords <- igraph::layout_with_kk(coords=theGrid, dim=somdim,
        igraph::mst(
          igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
            as.matrix(stats::dist(map$codes))^emcoords.pow)))
    } else if(emcoords=='fsom-mst') {
      emcoords <- igraph::layout_with_kk(dim=somdim, igraph::mst(
          igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
            as.matrix(stats::dist(map$codes)))))
    } else if(emcoords=='tsne') {
      emcoords <- Rtsne::Rtsne(map$codes, dims=somdim)$Y
    } else if(emcoords=='uwot::umap') {
      emcoords <- uwot::umap(map$codes, n_components=somdim)
    } else if(emcoords=='umap') {
      cfg <- umap::umap.defaults
      cfg$n_components <- somdim
      emcoords <- umap::umap(map$codes, cfg)$layout
    } else stop("unsupported emcoords method")
  }

  if(somdim==2) {
    if(dim(emcoords)[1] != x*y || dim(emcoords)[2] != 2) {
      stop("Embedding coordinates need to be of dimension (xdim*ydim, 2).")
    }
  } else {
    if(dim(emcoords)[1] != x*y*z || dim(emcoords)[2] != 3) {
      stop("Embedding coordinates need to be of dimension (xdim*ydim*zdim, 3).")
    }
  }

  codes <- t(map$codes)
  emcoords <- t(emcoords)

  embedding <- matrix(0, nrow=nrow(data), ncol=somdim)

  res <- .C("C_embedSOM",
    psomdim=as.integer(somdim),
    pn=as.integer(ndata),
    pdim=as.integer(dim),

    pboost=as.single(boost),
    pneighbors=as.integer(k),
    padjust=as.single(adjust),
    
    # the function now relies on the grid being arranged
    # reasonably, indexed row-by-row. If that changes, it is
    # necessary to pass in the map$grid as well.
    pncodes=as.integer(ncodes),

    points=as.single(points),
    koho=as.single(codes),
    emcoords=as.single(emcoords),

    embedding=as.single(embedding))

  matrix(res$embedding, nrow=nrow(data), ncol=somdim)
}
