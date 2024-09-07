#' Project velocities onto an embedding
#'
#' Project the velocity vector for each cell onto an existing low-dimensional embedding.
#'
#' @param x A numeric matrix of low-dimensional coordinates, e.g., after t-SNE.
#' Alternatively, a \linkS4class{SingleCellExperiment} containing such coordinates in its \code{\link{reducedDims}}.
#' @param vobj A \linkS4class{SingleCellExperiment} containing the output of the velocity calculations,
#' typically after running \code{\link{scvelo}}.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the ANY method, further arguments to pass to the \code{velocity_embedding} Python function from \pkg{scVelo}.
#'
#' For the SingleCellExperiment method, further arguments to pass to the ANY method.
#' @param use.dimred String or integer scalar specifying the reduced dimensions to retrieve from \code{x}.
#'
#' @details
#' This is a simple wrapper around the \code{scvelo.tools.velocity_embedding} function.
#' Briefly, we construct a cell-cell transition matrix where a cell is more likely to transition to one of its neighbors
#' if its velocity vector is pointing in the same direction as that neighbor.
#' The resulting matrix is then used to compute a weighted average of the positions in \code{x},
#' allowing us to compute a velocity in the low-dimensional embedding.
#' 
#' @return
#' A numeric matrix of the same dimensions as \code{x}, containing the projected velocity vectors in that embedding.
#'
#' @author Aaron Lun
#' @examples
#' out <- run_scvelo_example()
#'
#' # Making up a new embedding.
#' tsne.results <- matrix(rnorm(2*ncol(out)), ncol=2)
#'
#' # Projecting the future state of each cell:
#' projected <- embedVelocity(tsne.results, out)
#'
#' @export
#' @name embedVelocity
NULL

#' @importFrom SingleCellExperiment reducedDim<-
.embed_velocity <- function(x, vobj, ...) {
    reducedDim(vobj, "X_target") <- as.matrix(x)
    colnames(reducedDim(vobj, "X_target")) <- NULL
    basiliskRun(env=velo.env, fun=.run_embedder, vobj=vobj, ...)
}

#' @importFrom reticulate import
#' @importFrom zellkonverter SCE2AnnData
.run_embedder <- function(vobj, ...) {
    scv <- import("scvelo")

    args <- list(..., basis="target", autoscale=FALSE)
    adata <- SCE2AnnData(vobj)

    do.call(scv$tl$velocity_embedding, c(list(data=adata), args))

    adata$obsm["velocity_target"]
}

#' @export
#' @rdname embedVelocity
setGeneric("embedVelocity", function(x, vobj, ...) standardGeneric("embedVelocity"))

#' @export
#' @rdname embedVelocity
setMethod("embedVelocity", "ANY", .embed_velocity)

#' @export
#' @rdname embedVelocity
setMethod("embedVelocity", "SingleCellExperiment", function(x, vobj, ..., use.dimred=1) {
    .embed_velocity(reducedDim(x, use.dimred), vobj, ...)
})
