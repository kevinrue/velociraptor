# scvelo 0.3.1 and latest version of dependencies in zellkonverter::AnnDataDependencies() succesfully tested.
.scvelo_dependencies <- c(
    "scvelo==0.3.1",
    "anndata==0.10.2",
    "h5py==3.10.0",
    "hdf5==1.14.2",
    "natsort==8.4.0",
    "numpy==1.26.0",
    "packaging==23.2",
    "pandas==2.1.1",
    "python==3.11.5",
    "scipy==1.11.3"
)

#' @importFrom basilisk BasiliskEnvironment
#' @importFrom zellkonverter AnnDataDependencies
velo.env <- BasiliskEnvironment("env", "velociraptor",
    packages=.scvelo_dependencies, channels = c("bioconda", "conda-forge"))
