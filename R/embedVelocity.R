#' Project velocities onto an embedding
#'
#' Project the velocity vector for each cell onto an existing low-dimensional embedding.
#'
#' @param x A numeric matrix of low-dimensional coordinates, e.g., after t-SNE.
#' @param v A \linkS4class{SingleCellExperiment} containing the output of the velocity calculations,
#' typically after running \code{\link{scvelo}}.
#' @param ... Further arguments to pass to the \code{velocity_embedding} Python function from \pkg{scVelo}.
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
#' example(scvelo, echo=FALSE) # recycling that example.
#'
#' # Making up a new embedding.
#' tsne.results <- matrix(rnorm(2*ncol(out)), ncol=2)
#'
#' # Projecting the future state of each cell:
#' projected <- embedVelocity(tsne.results, out)
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
embedVelocity <- function(x, v, ...) {
    reducedDim(v, "X_target") <- as.matrix(x)
    basiliskRun(env=velo.env, fun=.run_embedder, v=v, ...)
}

#' @importFrom reticulate import
#' @importFrom zellkonverter SCE2AnnData
.run_embedder <- function(v, ...) {
    scv <- import("scvelo")

    args <- list(..., basis="target")
    adata <- SCE2AnnData(v)

    do.call(scv$tl$velocity_embedding, c(list(data=adata), args))

    adata$obsm["velocity_target"]
}
