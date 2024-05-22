# scvelo 0.2.5 and latest version of dependencies in zellkonverter::AnnDataDependencies() succesfully tested.
.scvelo_dependencies <- c(
    "scvelo=0.2.5",
    "matplotlib=3.7.3", # https://stackoverflow.com/questions/77128061/ydata-profiling-profilereport-attributeerror-module-matplotlib-cbook-has-no
    "tqdm=4.66.4", # required for progress bar
    "ipywidgets=8.1.2", # required for progress bar
    "jupyterlab=4.2.0" # required for progress bar
)

#' @importFrom basilisk BasiliskEnvironment
#' @importFrom zellkonverter AnnDataDependencies
velo.env <- BasiliskEnvironment("env", "velociraptor",
    packages=.scvelo_dependencies, channels = c("conda-forge", "bioconda"))
