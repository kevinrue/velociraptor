# scvelo 0.2.2 and latest version of dependencies in zellkonverter::.AnnDataDependencies succesfully tested.
.scvelo_dependencies <- c(
    "scvelo==0.2.2",
    "anndata==0.7.4",
    "h5py==2.10.0",
    "hdf5==1.10.6",
    "natsort==7.0.1",
    "numpy==1.19.1",
    "packaging==20.4",
    "pandas==1.1.2",
    "scipy==1.5.2"
)

#' @importFrom basilisk BasiliskEnvironment
#' @importFrom zellkonverter .AnnDataDependencies
velo.env <- BasiliskEnvironment("env", "velociraptor",
    packages=.scvelo_dependencies, channels = c("bioconda", "conda-forge"))
