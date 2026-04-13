# Project velocities onto an embedding

Project the velocity vector for each cell onto an existing
low-dimensional embedding.

## Usage

``` r
embedVelocity(x, vobj, ...)

# S4 method for class 'ANY'
embedVelocity(x, vobj, ...)

# S4 method for class 'SingleCellExperiment'
embedVelocity(x, vobj, ..., use.dimred = 1)
```

## Arguments

- x:

  A numeric matrix of low-dimensional coordinates, e.g., after t-SNE.
  Alternatively, a SingleCellExperiment containing such coordinates in
  its `reducedDims`.

- vobj:

  A SingleCellExperiment containing the output of the velocity
  calculations, typically after running
  [`scvelo`](https://kevinrue.github.io/velociraptor/reference/scvelo.md).

- ...:

  For the generic, further arguments to pass to specific methods.

  For the ANY method, further arguments to pass to the
  `velocity_embedding` Python function from scVelo.

  For the SingleCellExperiment method, further arguments to pass to the
  ANY method.

- use.dimred:

  String or integer scalar specifying the reduced dimensions to retrieve
  from `x`.

## Value

A numeric matrix of the same dimensions as `x`, containing the projected
velocity vectors in that embedding.

## Details

This is a simple wrapper around the `scvelo.tools.velocity_embedding`
function. Briefly, we construct a cell-cell transition matrix where a
cell is more likely to transition to one of its neighbors if its
velocity vector is pointing in the same direction as that neighbor. The
resulting matrix is then used to compute a weighted average of the
positions in `x`, allowing us to compute a velocity in the
low-dimensional embedding.

## Author

Aaron Lun

## Examples

``` r
example(scvelo, echo=FALSE) # recycling that example.
#> Loading required package: SingleCellExperiment
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> For native R and reading and writing of H5AD files, an R <AnnData> object, and
#> conversion to <SingleCellExperiment> or <Seurat> objects, check out the
#> anndataR package:
#> ℹ Install it from Bioconductor with `BiocManager::install("anndataR")`
#> ℹ See more at <https://bioconductor.org/packages/anndataR/>
#> This message is displayed once per session.
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")

# Making up a new embedding.
tsne.results <- matrix(rnorm(2*ncol(out)), ncol=2)

# Projecting the future state of each cell:
projected <- embedVelocity(tsne.results, out)
#> ℹ Using the 'X' assay as the X matrix
```
