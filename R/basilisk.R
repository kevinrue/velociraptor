#' @importFrom basilisk BasiliskEnvironment
#' @importFrom zellkonverter .AnnDataDependencies
velo.env <- BasiliskEnvironment("env", "velociraptor",
    packages=zellkonverter::.AnnDataDependencies, pip="scvelo==0.2.1")
