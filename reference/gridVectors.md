# Summarize vectors into a grid

Summarize the velocity vectors into a grid, usually for easy plotting.

## Usage

``` r
gridVectors(x, embedded, ...)

# S4 method for class 'ANY'
gridVectors(
  x,
  embedded,
  resolution = 40,
  scale = TRUE,
  as.data.frame = TRUE,
  return.intermediates = FALSE
)

# S4 method for class 'SingleCellExperiment'
gridVectors(x, embedded, ..., use.dimred = 1)
```

## Arguments

- x:

  A numeric matrix of low-dimensional coordinates, e.g., after t-SNE.
  Alternatively, a SingleCellExperiment containing such coordinates in
  its `reducedDims`.

- embedded:

  A low-dimensional projection of the velocity vectors into the
  embedding of `x`. This should be of the same dimensions as `x` and is
  typically produced by
  [`embedVelocity`](https://kevinrue.github.io/velociraptor/reference/embedVelocity.md).

- ...:

  For the generic, further arguments to pass to specific methods.

  For the SingleCellExperiment method, further arguments to pass to the
  ANY method.

- resolution:

  Integer scalar specifying the resolution of the grid, in terms of the
  number of grid intervals along each axis.

- scale:

  Logical scalar indicating whether the averaged vectors should be
  scaled by the grid resolution.

- as.data.frame:

  Logical scalar indicating whether the output should be a data.frame.
  If `FALSE`, a list of two matrices is returned.

- return.intermediates:

  Logical scalar indicating whether intermediate objects should also be
  returned. This enforces `as.data.frame=FALSE` and throws a warning if
  it is `TRUE`.

- use.dimred:

  String or integer scalar specifying the reduced dimensions to retrieve
  from `x`.

## Value

If `as.data.frame=FALSE`, a list is returned containing `start` and
`end`, two numeric matrices with one row per non-empty block in the grid
and one column per column in `x`. `start` contains the mean location of
all cells inside that block, and `end` contains the endpoint after
adding the (scaled) average of the block's cell's velocity vectors.

If `as.data.frame=TRUE`, a data.frame is returned with numeric columns
of the same contents as the list above. Column names are prefixed by
`start.*` and `end.*`.

If `return.intermediates=TRUE`, a list is returned (irrespective of the
value of `as.data.frame`) that in addition to `start` and `end` also
contains intermediate objects `limits` (the ranges in x and y), `delta`
(the grid intervals in x and y), `categories` (a DataFrame with integer
row and column indices for each cell that specify the grid field that it
is contained in), `grp` (numerical index of grid fields for each cell)
and `vec` (velocity vectors for non-empty grid fields).

## Details

This partitions the bounding box of `x` into a grid with `resolution`
units in each dimension. The locations and vectors of all cells in each
block are averaged to obtain a representative of that block. This is
most obviously useful for visualization to avoid overplotting of
velocity vectors.

If `scale=TRUE`, per-block vectors are scaled so that the median vector
length is comparable to the spacing between blocks. This improves
visualization when the scales of `x` and `embedded` are not immediately
comparable.

## See also

[`embedVelocity`](https://kevinrue.github.io/velociraptor/reference/embedVelocity.md),
to generate `embedded`.

## Author

Aaron Lun

## Examples

``` r
tsne.results <- matrix(rnorm(10000), ncol=2)
tsne.vectors <- matrix(rnorm(10000), ncol=2)

out <- gridVectors(tsne.results, tsne.vectors)

# Demonstration for plotting.
plot(tsne.results[,1], tsne.results[,2], col='grey')
arrows(out$start.1, out$start.2, out$end.1, out$end.2, length=0.05)

```
