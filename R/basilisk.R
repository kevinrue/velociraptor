#' @importFrom basilisk BasiliskEnvironment
velo.env <- BasiliskEnvironment("env", "velociraptor",
    packages=character(0), pip="scvelo==0.2.1")

#' @importFrom basilisk BasiliskEnvironment
velocyto.env <- BasiliskEnvironment("velocyto", "velociraptor",
    packages=c("numpy==1.18.5", "cython==0.29.20", "scipy==1.4.1", 
        "numba==0.49.1", "matplotlib==3.2.2", "scikit-learn==0.23.1",
        "h5py==2.10.0", "click==7.1.2"),
    pip="velocyto==0.17.17")
