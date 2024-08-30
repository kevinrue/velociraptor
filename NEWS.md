# velociraptor 1.15.8

* Update Conda environment using micromamba for MacOSX Arm.
* Fix switch between MacOSX and MacOSX Arm environments.

# velociraptor 1.15.7

* Update Conda environment using micromamba for Windows.

# velociraptor 1.15.6

* Update Conda environment to use `anaconda` channel on Linux.
  Passed GitHub Action <https://github.com/kevinrue/velociraptor/actions/runs/10612115572/job/29413105915>.

# velociraptor 1.15.5

* Update Conda environment to `svelo==0.3.2` on Linux. 

# velociraptor 1.15.4

* Set scvelo version triggering deprecation error to `0.3.1`.

# velociraptor 1.15.3

* Revert environment for Linux to the one of Bioconductor release 3.18.

# velociraptor 1.15.2

* Add environment for macOS (Intel); same environment as macos (M1).

# velociraptor 1.15.1

* Fix issue #63.
* Update `scvelo` to 0.3.2 (conda-forge) for macOS (M1) and Linux.
* Update `scvelo` to 0.2.5 (bioconda) for Windows.
* Add mechanism to switch Conda environment (and scvelo version) based on operating system and architecture.
* Use `scanpy.pp.neighbors` to calculate neighbors due to deprecation of automatic neighbor calculation in `scvelo.pp.moments`.
* Update vignette to document the change of default value for `n_neighbors` from scvelo (30) to scanpy (15).

# velociraptor 1.13.1

* Robust fallback mechanism using `basiliskRun` option `testload=`.

# velociraptor 1.9.3

* Pin python version in conda environment.

# velociraptor 1.9.2

* Remove column names from reduced dimension matrix in `gridVectors()`.

# velociraptor 1.5.2

* Remove column names of reduced dimension representation before velocity embedding.

# velociraptor 1.5.1

* Add example for `scvelo.params` argument.

# velociraptor 1.3.1

* Add typing_extensions to environment.

# velociraptor 1.1.6

* Move sanity check vignette to `inst/`.

# velociraptor 1.1.5

* Add Michael Stadler to package authors.

# velociraptor 1.1.4

* Fix typo in documentation.

# velociraptor 1.1.3

* Add vignette subdirectory with sanity checks.

# velociraptor 1.1.2

* Add functions `plotVelocity` and `plotVelocityStream`.

# velociraptor 1.1.1

* Refresh cached environments.

# velociraptor 1.1.0

* Bioconductor release 1.1.0.

# velociraptor 0.99.9

* Converted various functions to S4 generics for easier use with `SingleCellExperiment` objects.

# velociraptor 0.99.8

* Trigger new build to repeat `ExperimentHub` download.

# velociraptor 0.99.7

* Delete empty line to force cache update. See <https://github.com/rubocop-hq/rubocop/pull/4342#issuecomment-305449759>.

# velociraptor 0.99.6

* Set `autoscale=FALSE` in the call to `scvelo` function `velocity_embedding` to avoid issue related to Qt and plotting.

# velociraptor 0.99.5

* Trigger new build to check if Windows issue resolved itself.

# velociraptor 0.99.4

* Trigger new build to check whether TIMEOUT issue on Windows is reproducible.

# velociraptor 0.99.3

* Explicitly declare all Conda dependencies for `scvelo`.

# velociraptor 0.99.2

* Add hexsticker.

# velociraptor 0.99.1

* Remove .Rproj file from git repository.

# velociraptor 0.99.0

* First submission to Bioconductor.
