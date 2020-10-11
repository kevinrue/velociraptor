#' Summarize vectors into a grid
#'
#' Summarize the velocity vectors into a grid, usually for easy plotting.
#'
#' @param x A numeric matrix of low-dimensional coordinates, e.g., after t-SNE.
#' Alternatively, a \linkS4class{SingleCellExperiment} containing such coordinates in its \code{\link{reducedDims}}.
#' @param embedded A low-dimensional projection of the velocity vectors into the embedding of \code{x}.
#' This should be of the same dimensions as \code{x} and is typically produced by \code{\link{embedVelocity}}.
#' @param resolution Integer scalar specifying the resolution of the grid,
#' in terms of the number of grid intervals along each axis.
#' @param scale Logical scalar indicating whether the averaged vectors should be scaled by the grid resolution.
#' @param as.data.frame Logical scalar indicating whether the output should be a data.frame.
#' If \code{FALSE}, a list of two matrices is returned.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SingleCellExperiment method, further arguments to pass to the ANY method.
#' @param use.dimred String or integer scalar specifying the reduced dimensions to retrieve from \code{x}.
#'
#' @details
#' This partitions the bounding box of \code{x} into a grid with \code{resolution} units in each dimension.
#' The locations and vectors of all cells in each block are averaged to obtain a representative of that block.
#' This is most obviously useful for visualization to avoid overplotting of velocity vectors.
#'
#' If \code{scale=TRUE}, per-block vectors are scaled so that the median vector length is comparable to the spacing between blocks.
#' This improves visualization when the scales of \code{x} and \code{embedded} are not immediately comparable.
#'
#' @return
#' If \code{as.data.frame=FALSE}, a list is returned containing \code{start} and \code{end},
#' two numeric matrices with one row per non-empty block in the grid and one column per column in \code{x}.
#' \code{start} contains the mean location of all cells inside that block,
#' and \code{end} contains the endpoint after adding the (scaled) average of the block's cell's velocity vectors.
#'
#' If \code{as.data.frame=TRUE}, a data.frame is returned with numeric columns of the same contents as the list above.
#' Column names are prefixed by \code{start.*} and \code{end.*}.
#'
#' @author Aaron Lun
#' @examples
#' tsne.results <- matrix(rnorm(10000), ncol=2)
#' tsne.vectors <- matrix(rnorm(10000), ncol=2)
#'
#' out <- gridVectors(tsne.results, tsne.vectors)
#'
#' # Demonstration for plotting.
#' plot(tsne.results[,1], tsne.results[,2], col='grey')
#' arrows(out$start.1, out$start.2, out$end.1, out$end.2, length=0.05)
#'
#' @seealso 
#' \code{\link{embedVelocity}}, to generate \code{embedded}.
#' @name gridVectors
NULL

#' @export
#' @importFrom S4Vectors selfmatch DataFrame
#' @importFrom stats median
.grid_vectors <- function(x, embedded, resolution=40, scale=TRUE, as.data.frame=TRUE) {
    limits <- apply(x, 2, range)
    intercept <- limits[1,]
    delta <- (limits[2,] - intercept)/resolution

    categories <- t((t(x) - intercept)/delta)
    storage.mode(categories) <- "integer"
    categories <- DataFrame(categories)
    grp <- selfmatch(categories)

    tab <- as.integer(table(grp))
    pos <- rowsum(x, grp)/tab
    vec <- rowsum(embedded, grp)/tab

    if (scale) {
        d <- sqrt(rowSums(vec^2))
        target <- sqrt(sum(delta^2))
        vec <- vec * target/median(d)/2 # divide by 2 to avoid running over into another block.
    }

    FUN <- if (as.data.frame) data.frame else list
    FUN(start=pos, end=pos + vec)
}

#' @export
#' @rdname gridVectors
setGeneric("gridVectors", function(x, embedded, ...) standardGeneric("gridVectors"))

#' @export
#' @rdname gridVectors
setMethod("gridVectors", "ANY", .grid_vectors) 

#' @export
#' @rdname gridVectors
setMethod("gridVectors", "SingleCellExperiment", function(x, embedded, ..., use.dimred=1) {
    .grid_vectors(reducedDim(x, use.dimred), embedded, ...)
})
