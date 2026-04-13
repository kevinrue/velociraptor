# RNA velocity with scVelo

Perform RNA velocity calculations with the scVelo package.

## Usage

``` r
scvelo(x, ...)

# S4 method for class 'ANY'
scvelo(
  x,
  subset.row = NULL,
  sf.X = NULL,
  sf.spliced = NULL,
  sf.unspliced = NULL,
  use.theirs = FALSE,
  mode = c("steady_state", "deterministic", "stochastic", "dynamical"),
  scvelo.params = list(),
  dimred = NULL,
  ncomponents = 30,
  BPPARAM = SerialParam(),
  BSPARAM = bsparam()
)

# S4 method for class 'SummarizedExperiment'
scvelo(
  x,
  ...,
  assay.X = "counts",
  assay.spliced = "spliced",
  assay.unspliced = "unspliced"
)

# S4 method for class 'SingleCellExperiment'
scvelo(x, ..., sf.X = sizeFactors(x), dimred = NULL, use.dimred = NULL)
```

## Arguments

- x:

  A named list of three matrices of the same dimensions where genes are
  in rows and cells are in columns. The list should contain `"spliced"`
  and `"unspliced"` entries containing spliced and unspliced counts,
  respectively. It should also contain an `"X"` entry containing the
  “usual” count matrix, see details below.

  Alternatively, a SummarizedExperiment object containing three such
  matrices among its assays.

- ...:

  For the generic, further arguments to pass to specific methods. For
  the SummarizedExperiment and SingleCellExperiment methods, further
  arguments to pass to the ANY method.

- subset.row:

  A character, integer or logical vector specifying the genes to use for
  the velocity calculations. Defaults to all genes.

- sf.X:

  A numeric vector containing size factors for usual count matrix.
  Defaults to
  [`librarySizeFactors`](https://rdrr.io/pkg/scuttle/man/librarySizeFactors.html)
  on the `"X"` matrix in `x`.

- sf.spliced:

  A numeric vector containing size factors for the spliced counts for
  each cell. Defaults to
  [`librarySizeFactors`](https://rdrr.io/pkg/scuttle/man/librarySizeFactors.html)
  on the `"spliced"` matrix in `x`.

- sf.unspliced:

  A numeric vector containing size factors for the unspliced counts for
  each cell. Defaults to
  [`librarySizeFactors`](https://rdrr.io/pkg/scuttle/man/librarySizeFactors.html)
  on the `"unspliced"` matrix in `x`.

- use.theirs:

  Logical scalar indicating whether scVelo's gene filtering and
  normalization should be used.

- mode:

  String specifying the method to use to estimate the transcriptional
  dynamics.

- scvelo.params:

  List of lists containing arguments for individual scVelo functions,
  see details below.

- dimred:

  A low-dimensional representation of the cells with number of rows
  equal to the number of cells in `x`, used to find the nearest
  neighbors.

- ncomponents:

  Numeric scalar indicating the number of principal components to
  obtain. Only used if `use.theirs=FALSE` and `dimred=NULL`.

- BPPARAM:

  A BiocParallelParam object specifying whether the PCA calculations
  should be parallelized. Only used if `use.theirs=FALSE` and
  `dimred=NULL`.

- BSPARAM:

  A BiocSingularParam object specifying which algorithm should be used
  to perform the PCA. Only used if `use.theirs=FALSE` and `dimred=NULL`.

- assay.X:

  An integer scalar or string specifying the assay of `x` containing the
  usual count matrix.

- assay.spliced:

  An integer scalar or string specifying the assay of `x` containing the
  spliced counts.

- assay.unspliced:

  An integer scalar or string specifying the assay of `x` containing the
  unspliced counts.

- use.dimred:

  String naming the entry of
  [`reducedDims`](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)`(x)`
  to use for nearest neighbor calculations. Ignored if `dimred` is
  supplied.

## Value

A
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
is returned containing the output of the velocity calculations. Of
particular interest are:

- the `velocity_pseudotime` field in the `colData`, containing the
  velocity pseudotime for each cell.

- the `velocity` entry of the `assays`, containing the velocity vectors
  for each cell.

The output will always have number of columns equal to the number of
cells supplied in `x`, though the number of rows will depend on whether
any subsetting (if `subset.row` is supplied) or feature selection (if
`use.theirs=TRUE`) was performed.

## Details

This function uses the scVelo Python package
(<https://pypi.org/project/scvelo/>) for RNA velocity calculations. The
main difference from the original velocyto approach is that the
dynamical model of scVelo does not rely on the presence of observed
steady-state populations, which should improve the reliability of the
velocity calculations in general applications.

For consistency with other Bioconductor workflows, we perform as many
standard steps in R as we can before starting the velocity calculations
with scVelo. This involves:

1.  Size factor-based normalization with `sf.*` values and
    [`normalizeCounts`](https://rdrr.io/pkg/scuttle/man/normalizeCounts.html).
    For `"X"`, log-transformation is performed as well, while for the
    others, only scaling normalization is performed.

2.  Subsetting all matrices to `subset.row`, most typically to a subset
    of interest, e.g., highly variable genes. Note that, if set, any
    subsetting is done *after* normalization so that library sizes are
    correctly computed.

3.  If `dimred=NULL`, the PCA step on the log-expression values derived
    from the `"X"` matrix, using the specified `BSPARAM` to obtain the
    first `ncomponents` PCs.

This allows us to guarantee that, for example, the log-expression matrix
of HVGs or the PCA coordinates are the same as that used in other
applications like clustering or trajectory reconstruction.

Nonetheless, one can set `use.theirs=TRUE` to directly use the entire
scVelo normalization and filtering pipeline. This ignores all of the
size factors arguments (`sf.*`), all of the PCA-related arguments
(`ncomponents`, `BSPARAM`) and `subset.row`. However, if a
low-dimensionality result is supplied via `dimred` or `use.dimred`, the
scVelo PCA will always be omitted.

Upon first use, this function will instantiate a conda environment
containing the scVelo package. This is done via the basilisk package -
see the documentation for that package for trouble-shooting.

## Comments on the three matrices

Strictly speaking, only the spliced and unspliced matrices are necessary
for the velocity calculations. However, it is often the case that the
spliced matrix is not actually the same as a “usual” count matrix (e.g.,
generated by summing counts across all exons for all mapped genes). This
is due to differences in the handling of ambiguous reads that map across
exon-intron boundaries, or to genomic regions that can be either exonic
or intronic depending on the isoform; the spliced count matrix is more
likely to exclude such reads.

We request the usual count matrix as the `"X"` entry of `x` to ensure
that the PCA and nearest neighbor detection in scVelo are done on the
same data as that used in other steps of the large analysis (e.g.,
clustering, visualization, trajectory reconstruction). In practice, if
the usual count matrix is not available, one can often achieve
satisfactory results by simply re-using the spliced count matrix as both
the `"X"` and `"spliced"` entries of `x`.

Note that if reduced dimensions are supplied in `dimred`, any `"X"`
entry is only used to create the AnnData object and is not used in any
actual calculations.

## Additional arguments to Python

Additional arguments to scVelo functions are provided via
`scvelo.params`. This is a named list where each entry is named after a
function and is itself a named list of arguments for that function. The
following function names are currently recognized:

- `"filter_and_normalize"`, for gene selection and normalization. This
  is not used unless `use.theirs=TRUE`.

- `"moments"`, for PCA and nearest neighbor detection. The PCA is not
  performed if `dimred` or `use.dimred` is already supplied.

- `"recover_dynamics"`

- `"velocity"`

- `"velocity_graph"`

- `"velocity_pseudotime"`

- `"latent_time"`

- `"velocity_confidence"`

See the scVelo documentation for more details about the available
arguments and the examples below for a syntax example.

## Supported operating systems and architectures

scVelo dependencies are pinned in a Conda environment to ensure
reproducibility.

Differences in packages and versions available from Conda require
different environments for different operating systems and
architectures. basilisk.utils is used to determine the operating system
and architecture of the computer used to run `scvelo()`, using to the
appropriate Conda environment.

As of the latest velociraptor update (24 May 2024):

- All environments:

  tqdm and ipywidgets are installed to suppress the message "Unable to
  create progress bar".

- Linux:

  scVelo v0.3.2 from conda-forge is used. This is the latest version
  available to date. libtiff is pinned to v4.5.1 and pillow is pinned to
  v10.0.0
  (<https://github.com/conda-forge/libtiff-feedstock/issues/104#issuecomment-2375893029>),
  scipy is pinned to v1.13.1
  (<https://github.com/theislab/scvelo/issues/1260>).

- Linux AArch64:

  scVelo v0.3.2 from conda-forge is used. This is the latest version
  available to date. libtiff is pinned to v4.5.1 and pillow is pinned to
  v10.0.0
  (<https://github.com/conda-forge/libtiff-feedstock/issues/104#issuecomment-2375893029>),
  scipy is pinned to v1.13.1
  (<https://github.com/theislab/scvelo/issues/1260>). Please note that
  this environment has not been validated yet; it is derived from the
  environment from Linux (above) and requires additional testing to
  identify a working environment before pinning all the packages and
  versions in the environment.

- macOS:

  scVelo v0.3.2 from conda-forge is used. This is the latest version
  available to date. scipy is pinned to v1.13.1
  (<https://github.com/theislab/scvelo/issues/1260>).

- macOS ARM:

  scVelo v0.3.2 from conda-forge is used. This is the latest version
  available to date. scipy is pinned to v1.13.1
  (<https://github.com/theislab/scvelo/issues/1260>).

- Windows:

  scVelo v0.2.5 from bioconda is used. Later versions of scVelo depend
  on jaxlib which is not supported on Windows
  (<https://github.com/google/jax/issues/438>). matplotlib is pinned to
  v3.6.3 (<https://github.com/scverse/scanpy/issues/2411>), pandas is
  pinned to v1.5.2
  (<https://stackoverflow.com/questions/76234312/importerror-cannot-import-name-is-categorical-from-pandas-api-types>),
  and numpy is pinned to v1.21.1
  (<https://github.com/theislab/scvelo/issues/1109>).

## References

Bergen, V., Lange, M., Peidli, S. et al. Generalizing RNA velocity to
transient cell states through dynamical modeling. Nat Biotechnol 38,
1408–1414 (2020). <https://doi.org/10.1038/s41587-020-0591-3>

## Author

Aaron Lun, Charlotte Soneson, Kevin Rue-Albrecht

## Examples

``` r
# Using mock data to demonstrate the process:
library(scuttle)
sce1 <- mockSCE()
sce2 <- mockSCE()

spliced <- counts(sce1)
unspliced <- counts(sce2)

out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced))
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

# make scvelo use 10 rather than the default 30 neighbors to compute moments for velocity estimation:
out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced),
              scvelo.params=list(neighbors=list(n_neighbors=10L)))
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
```
