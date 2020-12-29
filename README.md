
# EmbedSOM

[![CRAN status](https://www.r-pkg.org/badges/version/EmbedSOM)](https://cran.r-project.org/package=EmbedSOM)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/EmbedSOM)](https://cran.r-project.org/package=EmbedSOM)

Fast embedding and visualization of multidimensional datasets, originally intended for flow/mass cytometry data. Compatible with FlowSOM (https://github.com/SofieVG/FlowSOM).

You may read more about EmbedSOM in the research paper here:

> Miroslav Kratochvíl, Abhishek Koladiya, and Jiří Vondrášek. "[Generalized EmbedSOM on quadtree-structured self-organizing maps](https://f1000research.com/articles/8-2120)" *F1000Research* **8** (2019). [doi:10.12688/f1000research.21642.2](https://doi.org/10.12688/f1000research.21642.2)

## Installation of R module

Use `devtools`:

```r
devtools::install_github('exaexa/EmbedSOM')
```

## Usage

EmbedSOM works by aligning the points to a precomputed self-organizing map (SOM). The main function `EmbedSOM` takes the SOM and data, and returns a matrix with 2D point coordinates on each row.

Quick way to get a visualization of multidimensional points saved in rows of `d`:

```r
map <- EmbedSOM::SOM(d, xdim=20, ydim=20)  # compute the self-organizing map
e <- EmbedSOM::EmbedSOM(d, map) # compute 2D coordinates of points
EmbedSOM::PlotEmbed(e) # plot the result with density
```

There are some parameters that affect speed, precision and shape of the embedding. Use `?EmbedSOM::EmbedSOM` to explore them in the documentation.

## HOW-TOs

To get started quickly, you can have a look at the vignettes:

- [Basic embedding on a toy dataset](https://bioinfo.uochb.cas.cz/embedsom/vignettes/basic.html)
- [Visualization of single-cell cytometry data from A FCS file](https://bioinfo.uochb.cas.cz/embedsom/vignettes/landmarks.html)
- [Advanced visualization of cytometry data with pseudotime](https://bioinfo.uochb.cas.cz/embedsom/vignettes/time.html)
- [Embedding 3D animal skeleton pointclouds to 2D](https://bioinfo.uochb.cas.cz/embedsom/vignettes/bones.html)

#### How to save an embedding to an FCS file?

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
EmbedSOM::PlotEmbed(e1, fsom=fs1)
# repeat as needed for file2.fcs, etc.
```

#### What are the color parameters of PlotEmbed?

See documentation in `?PlotEmbed`. By default, `PlotEmbed` plots a simple colored representation of point density. If supplied with a FCS column name (or number), it uses the a color scale similar to ColorBrewer's RdYlBu (with improvements for transparent stuff) to plot a single marker expression. Parameters `red`, `green` and `blue` can be used to set column names (or numbers) to mix RGB color from marker expressions.

`PlotEmbed` optionally accepts parameter `col` with a vector of R colors, which, if provided, is just forwarded to the internal `plot` function. For example, use `col=rgb(0,0,0,0.2)` for transparent black points.

**New!** if you need to mix more nicer colors than just the default RGB, use `ExprColors`.

#### How to plot the gazillions of the tiny points faster?
#### How to reduce size (and loading time) of the PDFs that contain scatterplots?

Use scattermore: https://github.com/exaexa/scattermore

#### How to select point subsets from the embedding?

A pretty fast (and still precise) way to dissect the dataset is to run a metaclustering on SOM clusters, and map the result to the individual points:

```r
clusters <- cutree(k=10, hclust(method='average', dist(map$codes)))[map$mapping[,1]]
```

After that, the metaclusters can be plotted in the embedding. Because the clustering is related to the small SOM-defined "pre-clusters" rather than the individual points, it is necessary to map the individual points to these clusters first:

```r
EmbedSOM::PlotEmbed(e, clust=clusters)
```

After you choose a metacluster in the embedding, use the color scale to find its number, then filter the points in your dataset to the corresponding subset. This example selects the point subset in metacluster number `5` and `6`:

```r
d <- d[clusters %in% c(5,6), ]
```

#### How to produce and display a 3D embedding?

There is now support for 3D SOM grids and 3D embedding. You need the customized SOM function from EmbedSOM:

```r
map <- EmbedSOM::SOM(someData, xdim=8, ydim=8, zdim=8)
e <- EmbedSOM::EmbedSOM(data=someData, map=map)
```

`PlotEmbed` and other functions do not work on 3D points in `embed`, but you may use other libraries to display the plots. For example, the `plot3D` library:

```r
plot3D::scatter3D(x=e[,1], y=e[,2], z=e[,3])
```

Interactive rotatable and zoomable plots are provided by the `rgl` library:

```r
rgl::points3d(x=e[,1], y=e[,2], z=e[,3])
```

#### What to do if embedding takes too long?

You may use parallelized versions of the algorithms. Several functions (`SOM`, `GQTSOM`, `EmbedSOM`) support setting `parallel=T`, which enables parallel processing; you may fine-tune the number of used CPUs by setting e.g. `threads=5`.

For SOM training, you need to explicitly switch to the parallelizable batch version, using `batch=T, parallel=T`.

#### How to activate the SIMD support? (i.e. how to get even more speed?)

Additionally, EmbedSOM has support for SIMD-assisted computation of both SOM and the embedding. If your CPU can work with SSE4 instructions (almost every `amd64` (a.k.a. `x64` a.k.a. `x86_64`) CPU built after around 2013 can do that), just tell R to compile your code with correct C++ flags, and SOM+EmbedSOM computation should get faster! (the usual speedup is at least around 3x, depending on the CPU and dataset shape)

First, add a correct line to the R `Makevars` config file:
```sh
 $ cat ~/.R/Makevars
CXXFLAGS += -O3 -march=native
```

After reinstalling EmbedSOM, SIMD code will be used by default. Note that only the functions from EmbedSOM are affected, i.e. you will need to use `EmbedSOM::SOM` instead of `FlowSOM::SOM` and `BuildSOM` to get this speedup.
