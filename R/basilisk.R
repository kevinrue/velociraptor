## Windows ----

.scvelo.dependencies.Windows <- list(
  packages = c(
    'scvelo==0.3.3'
  )
)

## LinuxAarch64 ----

.scvelo.dependencies.LinuxAarch64 <- list(
  channels = c("conda-forge"),
  packages = c(
    'ipywidgets==8.1.2',
    'libtiff==4.5.1',
    'pillow==10.0.0',
    'scipy==1.13.1',
    'scvelo==0.3.3',
    'tqdm==4.66.5'
  )
)

## Linux ----

.scvelo.dependencies.Linux <- list(
  channels = c("conda-forge"),
  packages = c(
    'scvelo==0.3.3'
  )
)

## MacOSXArm ----

.scvelo.dependencies.MacOSXArm <- list(
  packages = c(
    'scvelo==0.3.3'
  )
)

## MacOSX ----

.scvelo.dependencies.MacOSX <- list(
  channels = c("conda-forge", "bioconda"),
  packages = c(
    'scvelo==0.3.3'
  )
)

# Switch environment ----

#' @importFrom basilisk isWindows
#' @importFrom basilisk isLinuxAarch64
#' @importFrom basilisk isLinux
#' @importFrom basilisk isMacOSXArm
#' @importFrom basilisk isMacOSX
if (basilisk::isWindows()) {
  .scvelo_dependencies <- .scvelo.dependencies.Windows
} else if (basilisk::isLinuxAarch64()) {
  .scvelo_dependencies <- .scvelo.dependencies.LinuxAarch64
} else if (basilisk::isLinux()) {
  .scvelo_dependencies <- .scvelo.dependencies.Linux
} else if (basilisk::isMacOSXArm()) {
  .scvelo_dependencies <- .scvelo.dependencies.MacOSXArm
} else if (basilisk::isMacOSX()) {
  .scvelo_dependencies <- .scvelo.dependencies.MacOSX
} else {
  stop("Unsupported operating system or architecture.\n  Please open an issue at <https://github.com/kevinrue/velociraptor/issues> to request support.")
}

#' @importFrom basilisk BasiliskEnvironment
#' @importFrom zellkonverter AnnDataDependencies
velo.env <- BasiliskEnvironment("env", "velociraptor", packages=.scvelo_dependencies$packages)
