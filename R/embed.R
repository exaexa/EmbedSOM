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
#' @param emcoords Either a (xdim*ydim, 2)-matrix of 2D embedded coordinates of the SOM vertices, or one of 'flat' (default behavior), 'som' (adjust the SOM coords according to U-matrix distances), 'mst' (embed to MST-like structure), or 'fsom-mst' (embed to MST that should look exactly like that of FlowSOM)
#' @param emcoords.pow Exaggeration factor (power) of the distances in U-matrix used for some methods of auto-generating emcoords; default 1.
#' @return matrix with 2-D coordinates of the embedded cels
#'
#' @useDynLib EmbedSOM, .registration = TRUE
#' @export

EmbedSOM <- function(fsom=NULL, smooth=NULL, k=NULL, adjust=NULL,
                     data=NULL, map=NULL, importance=NULL,
                     emcoords='flat', emcoords.pow=1) {
  #TODO validate the sizes of data, colsUsed and codes.

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

  if(is.null(smooth))
    smooth <- 0

  if(is.null(k))
    k <- as.integer(1+sqrt(map$xdim*map$ydim))

  if(is.null(adjust)) {
    adjust <- 1
  }

  ndata <- nrow(data)
  colsUsed <- map$colsUsed
  if(is.null(colsUsed)) colsUsed <- (1:ncol(data))
  dim <- length(colsUsed)

  x <- map$xdim
  y <- map$ydim

  if (smooth< -30) {
    stop("Value of smooth must be at least -30.")
  }

  boost <- ((1+sqrt(5))/2)^(smooth-2)

  if (k<3) {
    stop("Use at least 3 neighbors for sane results!")
  }

  if (k>(x*y)) {
    stop("Too many neighbors!")
  }

  if(adjust<0) {
    stop("adjust must not be negative!");
  }

  if(!is.null(importance))
    points <- t(data[,colsUsed] * rep(importance, each=nrow(data)))
  else
    points <- t(data[,colsUsed])

  if(length(emcoords)==1) {
    if(emcoords=='flat') {
      emcoords <- as.matrix(expand.grid(1:x, 1:y))
    } else if(emcoords=='som') {
      emcoords <- igraph::layout_with_kk(coords=as.matrix(expand.grid(1:x, 1:y)),
        igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
          as.matrix(igraph::as_adjacency_matrix(igraph::make_lattice(c(x,y))))
          *
          as.matrix(dist(map$codes))^emcoords.pow))
    } else if(emcoords=='mst') {
      emcoords <- igraph::layout_with_kk(coords=as.matrix(expand.grid(1:x, 1:y)),
        igraph::mst(
          igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
            as.matrix(stats::dist(fsom$map$codes))^emcoords.pow)))
    } else if(emcoords=='fsom-mst') {
      emcoords <- igraph::layout_with_kk(igraph::mst(
          igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
            as.matrix(stats::dist(fsom$map$codes)))))
    } else stop("unsupported emcoords method")
  }

  if(dim(emcoords)[1] != x*y || dim(emcoords)[2] != 2) {
    stop("Embedding coordinates need to be of dimension (xdim*ydim, 2).")
  }

  codes <- t(map$codes)
  emcoords <- t(emcoords)

  embedding <- matrix(0, nrow=nrow(data), ncol=2)

  res <- .C("C_embedSOM",
    pn=as.integer(ndata),
    pdim=as.integer(dim),

    pboost=as.single(boost),
    pneighbors=as.integer(k),
    padjust=as.single(adjust),
    
    # the function now relies on the grid being arranged
    # reasonably, indexed row-by-row. If that changes, it is
    # necessary to pass in the map$grid as well.
    pxdim=as.integer(x),
    pydim=as.integer(y),

    points=as.single(points),
    koho=as.single(codes),
    emcoords=as.single(emcoords),

    embedding=as.single(embedding))

  matrix(res$embedding, nrow=nrow(data), ncol=2)
}
