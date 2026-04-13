# velociraptor

*velociraptor* provides an R toolkit for single-cell velocity
computation.

## Bioconductor release status

| Branch | R CMD check | Last updated |
|:--:|:--:|:--:|
| [*devel*](http://bioconductor.org/packages/devel/bioc/html/velociraptor.md) | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/velociraptor.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/velociraptor) | ![Latest commit](http://bioconductor.org/shields/lastcommit/devel/bioc/velociraptor.svg) |
| [*release*](http://bioconductor.org/packages/release/bioc/html/velociraptor.md) | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/velociraptor.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/velociraptor) | \![[Latest commit](http://bioconductor.org/shields/lastcommit/release/bioc/velociraptor.svg) |

## Installation

*velociraptor* can be easily installed from
[Bioconductor](https://bioconductor.org/packages/velociraptor/) using
[`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html):

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("velociraptor")
# or also...
BiocManager::install("velociraptor", dependencies = TRUE)
```

Setting `dependencies = TRUE` should ensure that all packages, including
the ones in the `Suggests:` field of the `DESCRIPTION` file, are
installed - this can be essential if you want to reproduce the code in
the vignette, for example.
