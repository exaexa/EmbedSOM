
# EmbedSOM

Fast embedding for flow/mass cytometry data. Best used with FlowSOM (https://github.com/SofieVG/FlowSOM).

## Installation of R module

Use `devtools`:

```r
devtools::install_github('exaexa/EmbedSOM')
```

## Usage

EmbedSOM works by aligning the cells to the FlowSOM-defined SOM (viewed as a smooth manifold). The main function `EmbedSOM` takes the SOM (present in the `$map` in FlowSOM objects) and returns a matrix with 2D cell coordinates on each row.

Quick way to get something out:

```r
fs <- FlowSOM::ReadInput('Levine_13dim_cleaned.fcs', scale=TRUE, transform=TRUE, toTransform=c(1:13))
fs <- FlowSOM::BuildSOM(fs, xdim=16, ydim=16, colsToUse=c(1:13))
e <- EmbedSOM::EmbedSOM(fs) # compute 2D coordinates of cells
par(mfrow=c(2,1))
EmbedSOM::PlotEmbed(e, fsom=fs)
```

(The FCS file can be downloaded from EmbedSOM website at http://bioinfo.uochb.cas.cz/embedsom/)

## EmbedSOM parameters

- `smooth`: Factor that affects how much difference between close and far nodes to create in neighborhood approximation. Increase this to produce "smoother" but possibly more convoluted embedding. (This functionality was originally done by the `n` parameter which selected the width of Gaussian for estimating the scores. The results were rounder but produced more artifacts and poor separation; it was replaced by `boost` parameter temporarily which did the same thing in a more complicated way.) Default value `0` is adjusted to be good for most datasets; lower values (`-1`, `-2`,...) produce sharper embeddings with less focus on projection and possibly more SOM-related artifacts. Higher values (`1`,`2`,...) produce smoother embeddings that better explains the large-scale structure of data, but some (possibly insignificant) small-scale details may get smoothed out. Sensible values are between `-10` and `10`.
- `k`: how many nearest SOM vertices to take into account at all (information from the `k+1`-th nearest SOM vertex is discarded). Performance depends quadratically on `k`. Increase to produce a more precise and smooth embedding. Setting between 10 and 50 is usually a good choice.
- `adjust`: Negative power factor for reducing the effect of non-local relevance measure on the outcome. Use `0` for plain projection; values above `1` usually push cells closer to respective SOM vertex positions.
- `fsom`: the FlowSOM object to embed
- `map`: optional map to use (e.g. if not present in the `fsom` object, or for embedding with different map)
- `data`: raw data matrix to be embedded (eg. if `fsom` object is not present). Must contain only the used columns, i.e. usually you want to use something like `data=myMatrix[,colsToUse]`)
- `importance`: same as for FlowSOM::BuildSOM. The `importance` passed to BuildSOM and EmbedSOM should be the same to prevent embedding artifacts (using different values breaks the k-NN calculation)
- `emcoords`: arbitraty set of map codes positions in the embedding. This may be generated using other more complicated dimensionality reduction algorithms without any loss of performance; EmbedSOM includes autogeneration of these positions (try e.g. parameter values `'som'`, `'mst'`, etc.).

## HOW-TOs

#### How to save an embedding to FCS file?

Use `flowCore` functionality to add any information to a FCS. The following template saves the scaled FlowSOM object data as-is, together with the embedding:

```r
fs <- FlowSOM::ReadInput('original.fcs', scale=T, ...)
fs <- FlowSOM::BuildSOM(fs, ...)
e <- EmbedSOM::EmbedSOM(fs, ...)
flowCore::write.FCS(new('flowFrame',
	exprs=as.matrix(data.frame(fs$data,
				   embedsom1=e[,1],
				   embedsom2=e[,2]))),
	'original_with_embedding.fcs')
```

See `flowCore` documentation for information about advanced FCS-writing functionality, e.g. for column descriptions.

#### How to align the population positions in embedding of multiple files?

Train a SOM on an aggregate file, and use it to embed the individual files. It is important to always apply the same scaling and transformations on all input files.

```r
fs <- FlowSOM::ReadInput(
	FlowSOM::AggregateFlowFrames(c('file1.fcs', 'file2.fcs', ...),
				     cTotal=100000),
	scale=T, transform=...)
n <- length(fs$scaled.scale)-2
map <- FlowSOM::SOM(fs)

fs1 <- FlowSOM::ReadInput('file1.fcs',
	scale=T,
	scaled.scale=fs$scaled.scale[1:n],
	scaled.center=fs$scaled.center[1:n],
	transform=...)
e1 <- EmbedSOM::EmbedSOM(fs1, map=map)
EmbedSOM::PlotEmbed(e1, fsom=fs1, xdim=10, ydim=10)
# repeat as needed for file2.fcs, etc.
```

#### What are the color parameters of PlotEmbed?

Please see documentation in `?PlotEmbed`. By default, `PlotEmbed` plots a simple colored representation of cell density. If supplied with a FCS column name (or number), it uses the `matlab2` color scale to plot a single marker expression. Parameters `red`, `green` and `blue` can be used to set column names (or numbers) to mix a RGB color from marker expressions.

`PlotEmbed` optionally accepts parameter `col` with a vector of R colors, which, if provided, is just forwarded to the internal `plot` function. For example, use `col=rgb(0,0,0,0.2)` for transparent black cells.

#### How to select cell subsets from the embedding?

First, a clustering method is needed to define the subsets. Supposing you already have a FlowSOM object `fs` with the SOM built, you can run e.g. the FlowSOM metaclustering to generate 10 clusters:

```r
clusters <- FlowSOM::metaClustering_consensus(fs$map$codes, k=10)
```

After that, the metaclusters can be plotted in the embedding. Because the clustering is related to the small FlowSOM "pre-clusters" rather than cells, it is also necessary to use the information from `fs$map$mapping` for getting the cluster information to single cell level:

```r
EmbedSOM::PlotEmbed(e, fsom=fs, clust=clusters[fs$map$mapping[,1]])
```

After you choose a metacluster in the embedding, use the color scale to find its number, then filter the cells in `fs` to the corresponding subset. This example selects the cell subset in metacluster number `5`:

```r
fs$data <- fs$data[clusters[fs$map$mapping[,1]]==5,]
```

Note that you must rebuild the SOM and re-embed the cells to work with the updated `fs` object.

#### How to produce and display a 3D embedding?

There is now support for 3D SOM grids and 3D embedding. You need the customized SOM function from EmbedSOM:

```r
map <- EmbedSOM::SOM(someData, xdim=8, ydim=8, zdim=8)
embed <- EmbedSOM::EmbedSOM(data=someData, map=map)
```

`PlotEmbed` and other functions do not work on 3D `embed` data, but you may use other libraries to display the plots. For example the `plot3D` library:

```r
plot3D::scatter3D(x=e[,1], y=e[,2], z=e[,3])
```

Interactive rotatable and zoomable plots can be viewed using the `rgl` library:

```r
rgl::points3d(x=e[,1], y=e[,2], z=e[,3])
```
