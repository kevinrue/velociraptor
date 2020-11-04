<img src="man/figures/logo.png" align="right" alt="logo.png" width="180" />

# velociraptor

<!-- badges: start -->
[![R build status](https://github.com/kevinrue/velociraptor/workflows/build_check_deploy/badge.svg)](https://github.com/kevinrue/velociraptor/actions)
[![Codecov.io coverage status](https://codecov.io/github/kevinrue/velociraptor/coverage.svg?branch=master)](https://codecov.io/github/kevinrue/velociraptor)
[![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/kevinrue/velociraptor)](https://hub.docker.com/r/kevinrue/velociraptor)
<!-- badges: end -->

_velociraptor_ provides an R toolkit for single-cell velocity computation.

# Installation

_velociraptor_ can be easily installed from [Bioconductor](https://bioconductor.org/packages/velociraptor/) using `BiocManager::install()`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("velociraptor")
# or also...
BiocManager::install("velociraptor", dependencies = TRUE)
```

Setting `dependencies = TRUE` should ensure that all packages, including the ones in the `Suggests:` field of the `DESCRIPTION` file, are installed - this can be essential if you want to reproduce the code in the vignette, for example.
