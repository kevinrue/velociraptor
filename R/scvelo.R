#' RNA velocity calculations with \pkg{scvelo}
#'
#' Perform RNA velocity calculations with the \pkg{scvelo} package.
#'
#' @param x A list of two matrices of the same dimensions where genes are in rows and cells are in columns.
#' The first matrix should contain spliced counts and the second matrix should contain unspliced counts.
#' If \code{norm=TRUE}, both matrices should contain normalized expression values instead.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing two such matrices in its assays.
#' @param ... For the generic, further arguments to pass to specific methods.
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.spliced An integer scalar or string specifying the assay of \code{x} containing the spliced counts.
#' @param assay.unspliced An integer scalar or string specifying the assay of \code{x} containing the unspliced counts.
#' @param norm Logical scalar indicating whether the matrices in \code{x} are already normalized.
#' In such cases, \code{sf.spliced} and \code{sf.unspliced} are ignored.
#' @param sf.spliced A numeric vector containing size factors for the spliced counts for each cell.
#' Defaults to \code{\link{librarySizeFactors}} on the spliced matrix.
#' @param sf.unspliced A numeric vector containing size factors for the unspliced counts for each cell.
#' Defaults to \code{\link{librarySizeFactors}} on the unspliced matrix.
#' @param subset.row A character, integer or logical vector specifying the genes to use for the velocity calculations.
#' Defaults to all genes but is most typically set to a subset of interesting genes, e.g., highly variable genes.
#' Note that, if set, any subsetting is done \emph{after} normalization so that library sizes are correctly computed.
#' @param use.theirs Logical scalar indicating whether \pkg{scvelo}'s gene filtering and normalization should be used,
#' in which case \code{sf.spliced}, \code{sf.unspliced} and \code{subset.row} are ignored.
#' @param mode String specifying the method to use to estimate the transcriptional dynamics.
#' @param scvelo.params List of lists, providing arguments for scVelo functions.
#' @param dimred A low-dimensional representation of the cells (N x k matrix, where N is the number of cells) that will be used to find the nearest neighbors for the moment calculations.
#'
#' @details
#' This function uses the \pkg{scvelo} package (\url{https://pypi.org/project/scvelo/}) to perform RNA velocity calculations.
#' The main difference from \code{\link{velocyto}} is that it does not rely on the presence of observed steady state populations,
#' which may improve the reliability of the velocity calculations in general applications.
#'
#' For consistency with other Bioconductor packages, we perform scaling normalization and subsetting
#' before starting the velocity calculations with \pkg{scvelo}.
#' This allows us to guarantee that, e.g., the spliced log-expression matrix of HVGs is the same as that used in other
#' applications like clustering or trajectory reconstruction.
#' Nonetheless, one can set \code{use.theirs=TRUE} to directly use the entire \pkg{scvelo} normalization and filtering pipeline.
#'
#' Upon first use, this function will instantiate a conda environment containing the \pkg{scvelo} package.
#'
#' @return
#' A \linkS4class{SummarizedExperiment} is returned containing the output of the velocity calculations.
#' TODO: SOME DOCUMENTATION REQUIRED.
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

#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle normalizeCounts
.scvelo <- function(x, subset.row=NULL, norm=FALSE,
    sf.spliced=NULL, sf.unspliced=NULL,
    use.theirs=FALSE,
    mode=c('steady_state', 'deterministic', 'stochastic', 'dynamical'),
    get.velocities=TRUE,
    scvelo.params=list(), dimred=NULL)
{
    spliced <- x[[1]]
    unspliced <- x[[2]]
    if (!identical(as.integer(dim(spliced)), as.integer(dim(unspliced)))) {
        stop("matrices in 'x' must have the same dimensions")
    }

    # Can't be bothered figuring this out.
    if (norm && use.theirs) {
        stop("need raw counts to use the 'scvelo' native pipeline")
    }

    if (!use.theirs) {
        if (!norm) {
            spliced <- normalizeCounts(spliced, sf.spliced, log=FALSE)
            unspliced <- normalizeCounts(unspliced, sf.unspliced, log=FALSE)
        }

        if (!is.null(subset.row)) {
            spliced <- spliced[subset.row,,drop=FALSE]
            unspliced <- unspliced[subset.row,,drop=FALSE]
        }

        # Find nearest neighbors
        if (length(dimred) > 0) {
            nn <- findKNN(dimred, k=29)$index
            nn <- cbind(seq_len(nrow(nn)), nn) - 1
            mode(nn) <- "integer"
            neighbors <- list(indices=nn, params=list(n_neighbors=30, n_pcs=ncol(dimred)))
        }
    } else {
        neighbors <- NULL
    }

    mode <- match.arg(mode)
    output <- basiliskRun(env=velo.env, fun=.run_scvelo,
        spliced=spliced, unspliced=unspliced,
        use.theirs=use.theirs, mode=mode,
        scvelo.params=scvelo.params,
        dimred=dimred, neighbors=neighbors)

    output
}

#' @importFrom reticulate import
#' @importFrom DelayedArray is_sparse t
.run_scvelo <- function(spliced, unspliced, use.theirs=FALSE, mode='dynamical', scvelo.params=list(), dimred=NULL, neighbors=NULL) {
    spliced <- t(.make_np_friendly(spliced))
    unspliced <- t(.make_np_friendly(unspliced))
    if (!is.null(dimred)) {
        dimred <- .make_np_friendly(dimred)
    }

    and <- import("anndata")
    scv <- import("scvelo")
    adata <- and$AnnData(spliced, layers=list(spliced=spliced, unspliced=unspliced))
    adata$obs_names <- rownames(spliced)
    adata$var_names <- colnames(spliced)
    if (!use.theirs && !is.null(dimred)) {
        adata$obsm <- list(X_pca = dimred)
    }

    if (use.theirs) {
        do.call(scv$pp$filter_and_normalize, c(list(data=adata),
                                               scvelo.params$filter_and_normalize))
    }

    do.call(scv$pp$moments, c(list(data=adata),
                              scvelo.params$moments))

    if (mode=="dynamical") {
        do.call(scv$tl$recover_dynamics, c(list(data=adata),
                                           scvelo.params$recover_dynamics))
    }

    scvelo.params$velocity$mode <- mode
    do.call(scv$tl$velocity, c(list(data=adata),
                               scvelo.params$velocity))

    do.call(scv$tl$velocity_graph, c(list(data=adata),
                                     scvelo.params$velocity_graph))

    do.call(scv$tl$velocity_pseudotime, c(list(adata=adata),
                                          scvelo.params$velocity_pseudotime))

    if (mode=="dynamical") {
        do.call(scv$tl$latent_time, c(list(data=adata),
                                      scvelo.params$latent_time))
    }

    do.call(scv$tl$velocity_confidence, c(list(data=adata),
                                          scvelo.params$velocity_confidence))

    .scvelo_anndata2sce(adata)
}

#' @importFrom S4Vectors make_zero_col_DFrame
.scvelo_anndata2sce <- function(adata) {
    assays <- .extractor_python_dict(adata$layers, c("Mu", "Ms", "velocity"), single=TRUE, transpose=TRUE)
    coldata <- .extractor_python_dict(adata$obs, names(adata$obs))
    rowdata <- .extractor_python_dict(adata$var, names(adata$var))
    metadata <- .extractor_python_dict(adata$uns, c('velocity_params', 'velocity_graph', 'velocity_graph_neg'), single=TRUE)

    if (!is.null(rowdata)) {
        rowdata <- do.call(DataFrame, rowdata)
    } else {
        rowdata <- make_zero_col_DFrame(nrow(assays[[1]]))
        rownames(rowdata) <- rownames(assays[[1]])
    }

    if (!is.null(coldata)) {
        coldata <- do.call(DataFrame, coldata)
    } else {
        coldata <- make_zero_col_DFrame(ncol(assays[[1]]))
        rownames(coldata) <- colnames(assays[[1]])
    }

    SummarizedExperiment(assays, colData=coldata, rowData=rowdata,
        metadata=metadata)
}

#' @export
#' @rdname scvelo
setGeneric("scvelo", function(x, ...) standardGeneric("scvelo"))

#' @export
#' @rdname scvelo
setMethod("scvelo", "ANY", .scvelo)

#' @export
#' @rdname scvelo
#' @importFrom BiocGenerics sizeFactors
setMethod("scvelo", "SummarizedExperiment", function(x, ...,
    assay.spliced="counts", assay.unspliced="unspliced",
    sf.spliced=sizeFactors(x), sf.unspliced=NULL, dimred="PCA")
{
    if (!is.null(dimred) && is(x, "SingleCellExperiment") && dimred %in% reducedDimNames(x)) {
        dimred <- reducedDim(x, dimred)
    } else {
        dimred=NULL
    }

    scvelo(list(assay(x, assay.spliced), assay(x, assay.unspliced)), ..., sf.spliced=sf.spliced, sf.unspliced=sf.unspliced, dimred=dimred)
})
