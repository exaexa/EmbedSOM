
# EmbedSOM

Fast embedding for flow/mass cytometry data. Best used with FlowSOM (https://github.com/SofieVG/FlowSOM).

## Installation of R module

Use `devtools`:

	devtools::install_github('exaexa/EmbedSOM')

## Usage

EmbedSOM works by aligning the cells to the FlowSOM-defined SOM (viewed as a smooth manifold). The main function `EmbedSOM` takes the SOM (present in the `$map` in FlowSOM objects) and returns a matrix with 2D cell coordinates on each row.

Quick way to get something out:

	fs <- FlowSOM::ReadInput('Levine_13dim_cleaned.fcs', scale=TRUE, transform=TRUE, toTransform=c(1:13))
	fs <- FlowSOM::BuildSOM(fs, xdim=16, ydim=16, colsToUse=c(1:13))
	e <- EmbedSOM::EmbedSOM(fs) # compute 2D coordinates of cells
	par(mfrow=c(2,1))
	EmbedSOM::PlotEmbed(e, fsom=fs)

(The FCS file can be downloaded from EmbedSOM website at http://bioinfo.uochb.cas.cz/embedsom/)

## EmbedSOM parameters

- `n`: how many nearest SOM vertices of the cell are considered as significant for the approximation (distance of the `n`-th neighbor is taken as `sigma` of a normal distributon of a relevance measure of SOM neighbors)
- `k`: how many nearest SOM vertices to take into account at all (information from the `k+1`-th nearest SOM vertex is discarded)
- `a`: used as a negative power for reducing the effect of non-local relevance measure on the outcome
- `fsom`: the FlowSOM object to embed
- `map`: optional map to use (e.g. if not present in the `fsom` object, or for embedding with different map)
- `data`: raw data matrix to be embedded (eg. if `fsom` object is not present). Must contain only the used columns, i.e. usually you want to use something like `data=myMatrix[,colsToUse]`)
- `importance`: same as for FlowSOM::BuildSOM. The `importance` passed to BuildSOM and EmbedSOM should be the same to prevent embedding artifacts (using different values breaks the k-NN calculation)

