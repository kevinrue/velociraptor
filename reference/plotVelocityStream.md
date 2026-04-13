# Velocity stream plot in low-dimensional space

Plot velocities embedded into low-dimensional space as a stream plot.
Stream lines are lines that follow the gradient in the velocity field
and illustrate paths that cells could follow based on observed RNA
velocities.

## Usage

``` r
plotVelocityStream(
  sce,
  embedded,
  use.dimred = 1,
  color_by = "#444444",
  color.alpha = 0.2,
  grid.resolution = 60,
  scale = TRUE,
  stream.L = 10,
  stream.min.L = 0,
  stream.res = 4,
  stream.width = 8,
  color.streamlines = FALSE,
  color.streamlines.map = c("#440154", "#482576", "#414487", "#35608D", "#2A788E",
    "#21908C", "#22A884", "#43BF71", "#7AD151", "#BBDF27", "#FDE725"),
  arrow.angle = 8,
  arrow.length = 0.8
)
```

## Arguments

- sce:

  A
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing low-dimensional coordinates, e.g., after t-SNE, in
  its
  [`reducedDims`](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).

- embedded:

  A low-dimensional projection of the velocity vectors into the
  embedding of `sce`. This should be of the same dimensions as `sce` and
  is typically produced by
  [`embedVelocity`](https://kevinrue.github.io/velociraptor/reference/embedVelocity.md).

- use.dimred:

  String or integer scalar specifying the reduced dimensions to retrieve
  from `sce`.

- color_by:

  A character scalar specifying a column in `colData(sce)` to color
  cells in the phase graph. Alternatively, `color_by` can be set to a
  valid R color to be used to color cells.

- color.alpha:

  An integer scalar giving the transparency of colored cells. Possible
  values are between 0 (fully transparent) and 1.0 (opaque).

- grid.resolution:

  Integer scalar specifying the resolution of the grid, in terms of the
  number of grid intervals along each axis.

- scale:

  Logical scalar indicating whether the averaged vectors should be
  scaled by the grid resolution.

- stream.L:

  Integer scalar giving the typical length of a streamline
  low-dimensional space units.

- stream.min.L:

  A numeric scalar with the minimum length of segments to be shown.

- stream.res:

  Numeric scalar specifying the resolution of estimated streamlines
  (higher numbers increase smoothness of lines but also the time for
  computation).

- stream.width:

  A numeric scalar controlling the width of streamlines.

- color.streamlines:

  Logical scalar. If `TRUE` streamlines will be colored by local
  velocity. Arrows cannot be shown in that case.

- color.streamlines.map:

  A character vector specifying the color range used for mapping local
  velocities to streamline colors. The default is
  `viridisLite::viridis(11)`.

- arrow.angle, arrow.length:

  Numeric scalars giving the `angle` and `length` of arrowheads.

## Value

A `ggplot2` object with the streamline plot.

## Details

`grid.resolution` and `scale` are passed to
[`gridVectors`](https://kevinrue.github.io/velociraptor/reference/gridVectors.md),
which is used to summarized the velocity vectors into an initial grid. A
full regular grid is computed from that and used in
[`geom_streamline`](https://eliocamp.github.io/metR/reference/geom_streamline.html)
to calculate streamlines. The following arguments are passed to the
arguments given in parenthesis of
[`geom_streamline`](https://eliocamp.github.io/metR/reference/geom_streamline.html):
`stream.L` (`L`), `stream.res` (`res`), `stream.min.L` (`min.L`),
`arrow.angle` (`arrow.angle`) and `arrow.length` (`arrow.length`).
Streamlines are computed by simple integration with a forward Euler
method, and `stream.L` and `stream.res` are used to compute the number
of steps and the time interval between steps for the integration.
`stream.width` is multiplied with `..step..` estimated by
[`geom_streamline`](https://eliocamp.github.io/metR/reference/geom_streamline.html)
to control the width of streamlines.

## See also

[`gridVectors`](https://kevinrue.github.io/velociraptor/reference/gridVectors.md)
used to summarize velocity vectors into a grid (velocity field), the
ggplot2 package used for plotting,
[`geom_streamline`](https://eliocamp.github.io/metR/reference/geom_streamline.html)
in package metR used to calculate and add streamlines from the RNA
velocity field to the plot,
[`viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html)
for creation of color palettes.

## Author

Michael Stadler

## Examples

``` r
library(scuttle)
set.seed(42)
sce1 <- mockSCE(ncells = 100, ngenes = 500)
sce2 <- mockSCE(ncells = 100, ngenes = 500)

datlist <- list(X=counts(sce1), spliced=counts(sce1), unspliced=counts(sce2))

out <- scvelo(datlist, mode = "dynamical")
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

em <- embedVelocity(reducedDim(out, 1), out)[,1:2]
#> ℹ Using the 'X' assay as the X matrix

# https://github.com/kevinrue/velociraptor/issues/90

if (!basilisk::isMacOSX()){
  plotVelocityStream(out, em)
  plotVelocityStream(out, em, color.streamlines = TRUE)
}
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the velociraptor package.
#>   Please report the issue at <https://github.com/kevinrue/velociraptor/issues>.
#> Warning: The dot-dot notation (`..step..`) was deprecated in ggplot2 3.4.0.
#> ℹ Please use `after_stat(step)` instead.
#> ℹ The deprecated feature was likely used in the velociraptor package.
#>   Please report the issue at <https://github.com/kevinrue/velociraptor/issues>.

```
