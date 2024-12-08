
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

guessDim <- function(dim, map) {
  if(is.null(dim)) {
    if(is.null(map$grid)) dim <- 2
    else dim <- ncol(map$grid)
  }
  dim
}
  
guessDistMethod <- function(dist.method, map) {
  if(is.null(dist.method)) {
    dist.method <- 'euclidean'
    if(!is.null(map$distf)) {
      if(map$distf==1) dist.method <- 'manhattan'
      if(map$distf==3) dist.method <- 'chebyshev'
    }
  }
  dist.method
}

#' Create a map by randomly selecting points
#'
#' @param data Input data matrix, with individual data points in rows
#' @param k How many points to sample
#' @param coordsFn a function to generate embedding coordinates (default none)
#' @return map object (without the grid, if `coordsFn` was not specified)
#'
#' @examples
#' d <- iris[,1:4]
#' EmbedSOM::PlotEmbed(
#'   EmbedSOM::EmbedSOM(
#'     data = d,
#'     map = EmbedSOM::RandomMap(d, 30, EmbedSOM::GraphCoords())),
#'   pch=19, clust=iris[,5]
#' )
#' @export
RandomMap <- function(data, k, coordsFn) {
  map <- list(codes = data[sample(nrow(data),k),])
  if(missing(coordsFn)) map
  else coordsFn(map)
}

#' Create a map from k-Means clusters
#'
#' May give better results than 'RandomMap' on data where random sampling
#' is complicated.
#' This does not use actual kMeans clustering, but re-uses the batch version of
#' [SOM()] with tiny radius (which makes it work the same as kMeans). In
#' consequence, the speedup of SOM function is applied here as well. Additionally,
#' because we don't need that amount of clustering precision, parameters `batch=F, rlen=1'
#' may give a satisfactory result very quickly.
#'
#' @param data Input data matrix, with individual data points in rows
#' @param k How many points to sample
#' @param coordsFn a function to generate embedding coordinates (default none)
#' @param batch Use batch-SOM training (effectively kMeans, default TRUE)
#' @param ... Passed to [SOM()], useful e.g. for 'parallel=T' or 'rlen=5'
#' @return map object (without the grid, if coordsFn was not specified)
#'
#' @examples
#' d <- iris[,1:4]
#' EmbedSOM::PlotEmbed(
#'   EmbedSOM::EmbedSOM(
#'     data = d,
#'     map = EmbedSOM::kMeansMap(d, 10, EmbedSOM::GraphCoords())),
#'   pch=19, clust=iris[,5]
#' )
#' @export
kMeansMap <- function(data, k, coordsFn, batch=T, ...) {
  map <- SOM(data,
             xdim=k, ydim=1,
             negRadius=0, negAlpha=0,
             batch=batch, radiusA=c(0.001, 0.001),
             ...)
  map$xdim <- NULL # pretend there was no SOM
  map$ydim <- NULL
  map$grid <- NULL
  if(missing(coordsFn)) map
  else coordsFn(map)
}

#' Add tSNE-based coordinates to a map
#'
#' @param dim Dimension of the result (passed to `tSNEFn` as `dims`)
#' @param tSNEFn tSNE function to run (default [Rtsne::Rtsne])
#' @param ... passed to `tSNEFn`
#' @return a function that transforms the map, usable as `coordsFn` parameter
#' @export
tSNECoords <- function(dim=NULL, tSNEFn=Rtsne::Rtsne, ...) {
  function(map) {
    dim <- guessDim(dim, map)
    map$grid <- tSNEFn(map$codes, dims=dim, ...)$Y
    map
  }
}

#' Add UMAP-based coordinates to a map
#'
#' @param dim Dimension of the result (passed to `UMAPFn` as `n_components`)
#' @param UMAPFn UMAP function to run (default [umap::umap] configured by [umap::umap.defaults])
#' @return a function that transforms the map, usable as `coordsFn` parameter
#' @export
UMAPCoords <- function(dim=NULL, UMAPFn=NULL) {
  if(is.null(UMAPFn)) UMAPFn <- function(x, dim) {
    cfg <- umap::umap.defaults
    cfg$n_components = dim
    umap::umap(x, cfg)$layout
  }
  function(map) {
    dim <- guessDim(dim, map)
    map$grid <- UMAPFn(map$codes, dim)
    map
  }
}

#' Add UMAP-based coordinates to a map, using the 'uwot' package
#'
#' @param dim Dimension of the result (passed to `uwotFn` as `dims`)
#' @param uwotFn UMAP function to run (default [uwot::umap])
#' @param ... passed to `uwotFn`
#' @return a function that transforms the map, usable as `coordsFn` parameter
#' @export
uwotCoords <- function(dim=NULL, uwotFn=uwot::umap, ...) {
  function(map) {
    dim <- guessDim(dim, map)
    map$grid <- uwotFn(map$codes, n_components=dim, ...)
    map
  }
}

#' Add U-Matrix-optimized embedding coordinates to the map
#'
#' The map must already contain a SOM grid with corresponding `xdim`,`ydim` (possibly `zdim`)
#'
#' @param dim Dimension of the result (passed to `layoutFn`)
#' @param dist.method The method to compute distances, passed to [stats::dist()] as parameter `method`
#' @param distFn Custom transformation function of the distance matrix
#' @param layoutFn iGraph-compatible graph layouting function (default [igraph::layout_with_kk])
#' @return a function that transforms the map, usable as 'coordsFn' parameter
#' @export
UMatrixCoords <- function(dim=NULL, dist.method=NULL, distFn=function(x)x, layoutFn=igraph::layout_with_kk) {
  function(map) {
    if(is.null(map$grid) || is.null(map$xdim) || is.null(map$ydim))
      stop("UMatrixCoords requires a map created by standard SOM algorithm")

    dim <- guessDim(dim, map)
    dist.method <- guessDistMethod(dist.method, map)
    map$grid <- igraph::layout_with_kk(coords=as.matrix(map$grid), dim=dim,
      igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
        as.matrix(igraph::as_adjacency_matrix(
          igraph::make_lattice(
            if(dim==2) c(map$xdim,map$ydim)
            else c(map$xdim,map$ydim,map$zdim)
          )
        ))
        *
        distFn(as.matrix(stats::dist(map$codes, method=dist.method)))))
    map
  }
}

#' Add MST-style embedding coordinates to the map
#'
#' @param dim Dimension of the result (passed to layoutFn)
#' @param dist.method The method to compute distances, passed to [stats::dist()] as parameter `method`
#' @param distFn Custom transformation function of the distance matrix
#' @param layoutFn iGraph-compatible graph layouting function (default [igraph::layout_with_kk()])
#' @return a function that transforms the map, usable as `coordsFn` parameter
#' @export
MSTCoords <- function(dim=NULL, dist.method=NULL, distFn=function(x)x, layoutFn=igraph::layout_with_kk) {
  function(map) {
    dim <- guessDim(dim, map)
    dist.method <- guessDistMethod(dist.method, map)
    map$grid <- layoutFn(
      coords = if(is.null(map$grid)) NULL else as.matrix(map$grid),
      dim = dim,
      igraph::mst(
        igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
          distFn(as.matrix(stats::dist(map$codes, method=dist.method))))))
    map
  }
}

#' Add Kamada-Kawai-generated embedding coordinates to the map
#'
#' This uses a complete graph on the map codebook, which brings overcrowding
#' problems. It is therefore useful to transform the distances for avoiding that
#' (e.g. by exponentiating them slightly using distFn function).
#'
#' @param dim Dimension of the result (passed to `layoutFn`)
#' @param dist.method The method to compute distances, passed to [stats::dist()] as parameter `method`
#' @param distFn Custom transformation function of the distance matrix
#' @param layoutFn iGraph-compatible graph layouting function (default [igraph::layout_with_kk])
#' @return a function that transforms the map, usable as `coordsFn` parameter
#' @export
GraphCoords <- function(dim=NULL, dist.method=NULL, distFn=function(x)x, layoutFn=igraph::layout_with_kk) {
  function(map) {
    dim <- guessDim(dim, map)
    dist.method <- guessDistMethod(dist.method, map)
    map$grid <- layoutFn(
      dim=dim,
      igraph::graph_from_adjacency_matrix(mode='undirected', weighted=T,
        distFn(as.matrix(stats::dist(map$codes, method=dist.method)))))
    map
  }
}

#' Add KNN-topology-based embedding coordinates to the map
#'
#' Internally, this does not use [FNN::get.knn()] anymore.
#'
#' @param k Size of the neighborhoods (default 4)
#' @param dim Dimension of the result (passed to `layoutFn`)
#' @param dist.method The method to compute distances, passed to [stats::dist()] as parameter `method`
#' @param distFn Custom transformation function of the distance matrix
#' @param layoutFn iGraph-compatible graph layouting function (default [igraph::layout_with_kk])
#' @return a function that transforms the map, usable as `coordsFn` parameter
#' @export
kNNCoords <- function(k=4, dim=NULL, dist.method=NULL, distFn=function(x)x, layoutFn=igraph::layout_with_kk) {
  function(map) {
    dim <- guessDim(dim, map)
    dist.method <- guessDistMethod(dist.method, map)
    dist.matrix <- as.matrix(stats::dist(map$codes, method=dist.method))
    n <- nrow(map$codes)
    kns <- cbind(rep(1:n, each=k), as.vector(apply(dist.matrix, 1, function(x) order(x)[2:(k+1)])))
    adj <- matrix(0, n, n)
    adj[kns[,2]+n*(kns[,1]-1)] <- 1
    adj[kns[,1]+n*(kns[,2]-1)] <- 1
    adj <- adj*(distFn(dist.matrix))

    map$grid <- layoutFn(dim=dim,
      igraph::graph_from_adjacency_matrix(adj, mode='undirected', weighted=T))
    map
  }
}

