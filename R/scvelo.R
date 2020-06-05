#' RNA velocity calculations with \pkg{scvelo}
#'
#' Perform RNA velocity calculations with the \pkg{scvelo} package.
#'
#' @param x A list of two matrices of the same dimensions where genes are in rows and cells are in columns.
#' The first matrix should contain spliced counts and te second matrix should contain unspliced counts.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing two such matrices in its assays.
#' @param ... For the generic, further arguments to pass to specific methods.
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.spliced An integer scalar or string specifying the assay of \code{x} containing the spliced counts.
#' @param assay.unspliced An integer scalar or string specifying the assay of \code{x} containing the unspliced counts.
#'
#' @details
#' This function uses the \pkg{scvelo} package (\url{https://pypi.org/project/scvelo/}) to perform RNA velocity calculations.
#' The main difference from \code{\link{velocyto}} is that it does not rely on the presence of observed steady state populations.
#' 
#' Upon first use, this function will instantiate a conda environment containing the \pkg{scvelo} package.
#' 
#' @return 
#' By default, a \linkS4class{DataFrame} is returned containing:
#'
#' @examples
#' # Using mock data to demonstrate the process:
#' library(scuttle)
#' sce1 <- mockSCE()
#' sce2 <- mockSCE()
#'
#' spliced <- counts(sce1)
#' unspliced <- counts(sce2)
#'
#' out <- scvelo(list(spliced, unspliced))
#' 
#' @author Aaron Lun
#' @name scvelo
NULL

#' @importFrom reticulate import
.run_scvelo <- function(spliced, unspliced) {
    and <- import("anndata")
    scv <- import("scvelo")

    adata <- and$AnnData(spliced, layers=list(spliced=spliced, unspliced=unspliced))
    scv$pp$filter_and_normalize(adata)
    scv$pp$moments(adata) 

    scv$tl$velocity(adata, mode='stochastic')
    scv$tl$velocity_graph(adata)
    scv$tl$velocity_pseudotime(adata)

    list(pseudotime=adata$obs$velocity_pseudotime,
        roots=adata$obs$root_cells,
        endpoints=adata$obs$end_points)
}

#' @export
#' @rdname scvelo
setGeneric("scvelo", function(x, ...) standardGeneric("scvelo"))

#' @export
#' @rdname scvelo
#' @importFrom S4Vectors DataFrame
setMethod("scvelo", "ANY", function(x) {
    spliced <- x[[1]]
    unspliced <- x[[2]]
    if (!identical(as.integer(dim(spliced)), as.integer(dim(unspliced)))) {
        stop("matrices in 'x' must have the same dimensions")
    } 

    output <- basiliskRun(env=velo.env, fun=.run_scvelo, 
        spliced=spliced, unspliced=unspliced)

    do.call(DataFrame, output)
})

#' @export
#' @rdname scvelo
setMethod("scvelo", "SummarizedExperiment", function(x, ..., assay.spliced="spliced", assay.unspliced="unspliced") {
    x <- list(assay(x, assay.spliced), assay(x, assay.unspliced))
    scvelo(x, ...)
})
