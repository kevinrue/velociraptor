#' RNA velocity calculations with \pkg{scVelo}
#'
#' Perform RNA velocity calculations with the \pkg{scVelo} package.
#'
#' @param x A list of three matrices of the same dimensions where genes are in rows and cells are in columns.
#' The first matrix should contain the counts used for calculating reduced dimension representations and cell-to-cell distances.
#' The second matrix should contain spliced counts and the third matrix should contain unspliced counts.
#' If \code{norm=TRUE}, all matrices should contain normalized (but not log-transformed) expression values instead.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing three such matrices among its assays.
#' @param ... For the generic, further arguments to pass to specific methods.
#' For the SummarizedExperiment and SingleCellExperiment methods, further arguments to pass to the ANY method.
#' @param assay.x An integer scalar or string specifying the assay of \code{x} containing the counts that will be used for dimension reduction and distance calculations.
#' @param assay.spliced An integer scalar or string specifying the assay of \code{x} containing the spliced counts.
#' @param assay.unspliced An integer scalar or string specifying the assay of \code{x} containing the unspliced counts.
#' @param norm Logical scalar indicating whether the matrices in \code{x} are already normalized.
#' In such cases, \code{sf.x}, \code{sf.spliced} and \code{sf.unspliced} are ignored.
#' @param sf.x A numeric vector containing size factors for the counts in \code{x} for each cell.
#' Defaults to \code{\link{librarySizeFactors}} on the \code{x} count matrix.
#' @param sf.spliced A numeric vector containing size factors for the spliced counts for each cell.
#' Defaults to \code{\link{librarySizeFactors}} on the spliced matrix.
#' @param sf.unspliced A numeric vector containing size factors for the unspliced counts for each cell.
#' Defaults to \code{\link{librarySizeFactors}} on the unspliced matrix.
#' @param subset.row A character, integer or logical vector specifying the genes to use for the velocity calculations.
#' Defaults to all genes but is most typically set to a subset of interesting genes, e.g., highly variable genes.
#' Note that, if set, any subsetting is done \emph{after} normalization so that library sizes are correctly computed.
#' @param use.theirs Logical scalar indicating whether \pkg{scVelo}'s gene filtering and normalization should be used,
#' in which case \code{sf.counts}, \code{sf.spliced}, \code{sf.unspliced} and \code{subset.row} are ignored.
#' @param mode String specifying the method to use to estimate the transcriptional dynamics.
#' @param scvelo.params List of lists, providing arguments for scVelo functions.
#' @param dimred A low-dimensional representation of the cells (N x k matrix, where N is the number of cells) that will be used to find the nearest neighbors for the moment calculations.
#' Alternatively, if \code{x} is a \linkS4class{SingleCellExperiment}, a character scalar with the name of a reduced dimension representation in \code{x}.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain, if \code{use.theirs} is FALSE and no \code{dimred} is supplied.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA calculations should be parallelized.
#'
#' @details
#' This function uses the \pkg{scVelo} package (\url{https://pypi.org/project/scvelo/}) to perform RNA velocity calculations.
#' The main difference from \pkg{velocyto} is that the dynamical model of \pkg{scVelo} does not rely on the presence of observed steady state populations,
#' which may improve the reliability of the velocity calculations in general applications.
#'
#' For consistency with other Bioconductor packages, we perform scaling normalization and subsetting
#' before starting the velocity calculations with \pkg{scVelo}.
#' This allows us to guarantee that, e.g., the spliced log-expression matrix of HVGs is the same as that used in other
#' applications like clustering or trajectory reconstruction.
#' Nonetheless, one can set \code{use.theirs=TRUE} to directly use the entire \pkg{scVelo} normalization and filtering pipeline.
#'
#' Upon first use, this function will instantiate a conda environment containing the \pkg{scVelo} package.
#'
#' Additional arguments to \pkg{scVelo} functions are provided via the \code{scvelo.params} list.
#' The arguments for a given function are given as a named list.
#' The following function names are currently recognized:
#' \itemize{
#' \item filter_and_normalize
#' \item moments
#' \item recover_dynamics
#' \item velocity
#' \item velocity_graph
#' \item velocity_pseudotime
#' \item latent_time
#' \item velocity_confidence
#' }
#' For more information about the arguments of these functions, see the \pkg{scVelo} documentation.
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
#' out <- scvelo(list(spliced, spliced, unspliced), use.theirs=TRUE)
#'
#' @author Aaron Lun, Charlotte Soneson
#' @name scvelo
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle normalizeCounts
#' @importFrom BiocSingular bsparam
#' @importFrom BiocParallel SerialParam
.scvelo <- function(x, subset.row=NULL, norm=FALSE,
    sf.x=NULL, sf.spliced=NULL, sf.unspliced=NULL,
    use.theirs=FALSE,
    mode=c('steady_state', 'deterministic', 'stochastic', 'dynamical'),
    scvelo.params=list(), dimred=NULL, ncomponents=30, BPPARAM=SerialParam(),
    BSPARAM=bsparam())
{
    spliced <- x[[2]]
    unspliced <- x[[3]]
    x <- x[[1]]
    if (!all(identical(as.integer(dim(spliced)), as.integer(dim(unspliced))),
             identical(as.integer(dim(spliced)), as.integer(dim(x))))) {
        stop("matrices in 'x' must have the same dimensions")
    }

    # Can't be bothered figuring this out.
    if (norm && use.theirs) {
        stop("need raw counts to use the 'scvelo' native pipeline")
    }

    if (!use.theirs) {
        if (!norm) {
            x <- normalizeCounts(x, sf.x, log=TRUE)
            spliced <- normalizeCounts(spliced, sf.spliced, log=FALSE)
            unspliced <- normalizeCounts(unspliced, sf.unspliced, log=FALSE)
        } else {
            x <- log2(x + 1)
        }

        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
            spliced <- spliced[subset.row,,drop=FALSE]
            unspliced <- unspliced[subset.row,,drop=FALSE]
        }

        if (is.null(dimred)) {
            dimred <- BiocSingular::runPCA(t(x), rank=ncomponents,
                                           BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        }

    }

    mode <- match.arg(mode)
    output <- basiliskRun(env=velo.env, fun=.run_scvelo,
        x=x, spliced=spliced, unspliced=unspliced,
        use.theirs=use.theirs, mode=mode,
        scvelo.params=scvelo.params,
        dimred=dimred)

    output
}

#' @importFrom reticulate import
#' @importFrom DelayedArray is_sparse t
.run_scvelo <- function(x, spliced, unspliced, use.theirs=FALSE, mode='dynamical', scvelo.params=list(), dimred=NULL) {
    x <- t(.make_np_friendly(x))
    spliced <- t(.make_np_friendly(spliced))
    unspliced <- t(.make_np_friendly(unspliced))

    if (!is.null(dimred)) {
        dimred <- .make_np_friendly(dimred)
    }

    and <- import("anndata")
    scv <- import("scvelo")
    adata <- and$AnnData(x, layers=list(spliced=spliced, unspliced=unspliced))
    adata$obs_names <- rownames(spliced)
    adata$var_names <- colnames(spliced)

    ## A supplied dimred will be used even if use.theirs=TRUE
    if (!is.null(dimred)) {
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
    assay.x="counts", assay.spliced="spliced", assay.unspliced="unspliced",
    sf.x=sizeFactors(x), sf.spliced=NULL, sf.unspliced=NULL)
{
    scvelo(list(assay(x, assay.x), assay(x, assay.spliced), assay(x, assay.unspliced)), ..., sf.x=sf.x, sf.spliced=sf.spliced, sf.unspliced=sf.unspliced, dimred=NULL)
})

#' @export
#' @rdname scvelo
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
setMethod("scvelo", "SingleCellExperiment", function(x, ...,
    assay.x="counts", assay.spliced="spliced", assay.unspliced="unspliced",
    sf.x=sizeFactors(x), sf.spliced=NULL, sf.unspliced=NULL, dimred="PCA")
{
    if (!is.null(dimred) && dimred %in% reducedDimNames(x)) {
        dimred <- reducedDim(x, dimred)
    } else {
        dimred <- NULL
    }

    scvelo(list(assay(x, assay.x), assay(x, assay.spliced), assay(x, assay.unspliced)), ..., sf.x=sf.x, sf.spliced=sf.spliced, sf.unspliced=sf.unspliced, dimred=dimred)
})
