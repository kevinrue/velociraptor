#' RNA velocity with \pkg{scVelo}
#'
#' Perform RNA velocity calculations with the \pkg{scVelo} package.
#'
#' @param x A named list of three matrices of the same dimensions where genes are in rows and cells are in columns.
#' The list should contain \code{"spliced"} and \code{"unspliced"} entries containing spliced and unspliced counts, respectively.
#' It should also contain an \code{"X"} entry containing the \dQuote{usual} count matrix, see details below.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing three such matrices among its assays.
#' @param ... For the generic, further arguments to pass to specific methods.
#' For the SummarizedExperiment and SingleCellExperiment methods, further arguments to pass to the ANY method.
#' @param assay.X An integer scalar or string specifying the assay of \code{x} containing the usual count matrix.
#' @param assay.spliced An integer scalar or string specifying the assay of \code{x} containing the spliced counts.
#' @param assay.unspliced An integer scalar or string specifying the assay of \code{x} containing the unspliced counts.
#' @param sf.X A numeric vector containing size factors for usual count matrix.
#' Defaults to \code{\link{librarySizeFactors}} on the \code{"X"} matrix in \code{x}.
#' @param sf.spliced A numeric vector containing size factors for the spliced counts for each cell.
#' Defaults to \code{\link{librarySizeFactors}} on the \code{"spliced"} matrix in \code{x}.
#' @param sf.unspliced A numeric vector containing size factors for the unspliced counts for each cell.
#' Defaults to \code{\link{librarySizeFactors}} on the \code{"unspliced"} matrix in \code{x}.
#' @param subset.row A character, integer or logical vector specifying the genes to use for the velocity calculations.
#' Defaults to all genes.
#' @param use.theirs Logical scalar indicating whether \pkg{scVelo}'s gene filtering and normalization should be used.
#' @param mode String specifying the method to use to estimate the transcriptional dynamics.
#' @param scvelo.params List of lists containing arguments for individual \pkg{scVelo} functions, see details below.
#' @param dimred A low-dimensional representation of the cells with number of rows equal to the number of cells in \code{x},
#' used to find the nearest neighbors.
#' @param use.dimred String naming the entry of \code{\link{reducedDims}(x)} to use for nearest neighbor calculations.
#' Ignored if \code{dimred} is supplied.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' Only used if \code{use.theirs=FALSE} and \code{dimred=NULL}.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' Only used if \code{use.theirs=FALSE} and \code{dimred=NULL}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA calculations should be parallelized.
#' Only used if \code{use.theirs=FALSE} and \code{dimred=NULL}.
#'
#' @details
#' This function uses the \pkg{scVelo} Python package (\url{https://pypi.org/project/scvelo/}) for RNA velocity calculations.
#' The main difference from the original \pkg{velocyto} approach is that the dynamical model of \pkg{scVelo}
#' does not rely on the presence of observed steady-state populations,
#' which should improve the reliability of the velocity calculations in general applications.
#'
#' For consistency with other Bioconductor workflows, we perform as many standard steps in R as we can
#' before starting the velocity calculations with \pkg{scVelo}.
#' This involves:
#' \enumerate{
#' \item Size factor-based normalization with \code{sf.*} values and \code{\link{normalizeCounts}}.
#' For \code{"X"}, log-transformation is performed as well, while for the others, only scaling normalization is performed.
#' \item Subsetting all matrices to \code{subset.row}, most typically to a subset of interest, e.g., highly variable genes.
#' Note that, if set, any subsetting is done \emph{after} normalization so that library sizes are correctly computed.
#' \item If \code{dimred=NULL}, the PCA step on the log-expression values derived from the \code{"X"} matrix,
#' using the specified \code{BSPARAM} to obtain the first \code{ncomponents} PCs.
#' }
#' This allows us to guarantee that, for example, the log-expression matrix of HVGs or the PCA coordinates
#' are the same as that used in other applications like clustering or trajectory reconstruction.
#'
#' Nonetheless, one can set \code{use.theirs=TRUE} to directly use the entire \pkg{scVelo} normalization and filtering pipeline.
#' This ignores all of the size factors arguments (\code{sf.*}),
#' all of the PCA-related arguments (\code{ncomponents}, \code{BSPARAM}) and \code{subset.row}.
#' However, if a low-dimensionality result is supplied via \code{dimred} or \code{use.dimred},
#' the \pkg{scVelo} PCA will always be omitted.
#'
#' Upon first use, this function will instantiate a conda environment containing the \pkg{scVelo} package.
#' This is done via the \pkg{basilisk} package - see the documentation for that package for trouble-shooting.
#'
#' @section Comments on the three matrices:
#' Strictly speaking, only the spliced and unspliced matrices are necessary for the velocity calculations.
#' However, it is often the case that the spliced matrix is not actually the same as a \dQuote{usual} count matrix
#' (e.g., generated by summing counts across all exons for all mapped genes).
#' This is due to differences in the handling of ambiguous reads that map across exon-intron boundaries,
#' or to genomic regions that can be either exonic or intronic depending on the isoform;
#' the spliced count matrix is more likely to exclude such reads.
#'
#' We request the usual count matrix as the \code{"X"} entry of \code{x} to ensure that
#' the PCA and nearest neighbor detection in \pkg{scVelo} are done on the same data
#' as that used in other steps of the large analysis (e.g., clustering, visualization, trajectory reconstruction).
#' In practice, if the usual count matrix is not available, one can often achieve satisfactory results
#' by simply re-using the spliced count matrix as both the \code{"X"} and \code{"spliced"} entries of \code{x}.
#'
#' Note that if reduced dimensions are supplied in \code{dimred},
#' any \code{"X"} entry is only used to create the AnnData object and is not used in any actual calculations.
#'
#' @section Additional arguments to Python:
#' Additional arguments to \pkg{scVelo} functions are provided via \code{scvelo.params}.
#' This is a named list where each entry is named after a function and is itself a named list of arguments for that function.
#' The following function names are currently recognized:
#' \itemize{
#' \item \code{"filter_and_normalize"}, for gene selection and normalization.
#' This is not used unless \code{use.theirs=TRUE}.
#' \item \code{"moments"}, for PCA and nearest neighbor detection.
#' The PCA is not performed if \code{dimred} or \code{use.dimred} is already supplied.
#' \item \code{"recover_dynamics"}
#' \item \code{"velocity"}
#' \item \code{"velocity_graph"}
#' \item \code{"velocity_pseudotime"}
#' \item \code{"latent_time"}
#' \item \code{"velocity_confidence"}
#' }
#' See the \pkg{scVelo} documentation for more details about the available arguments and the examples below for a syntax example.
#'
#' @return
#' A \linkS4class{SingleCellExperiment} is returned containing the output of the velocity calculations.
#' Of particular interest are:
#' \itemize{
#' \item the \code{velocity_pseudotime} field in the \code{\link{colData}},
#' containing the velocity pseudotime for each cell.
#' \item the \code{velocity} entry of the \code{\link{assays}},
#' containing the velocity vectors for each cell.
#' }
#' The output will always have number of columns equal to the number of cells supplied in \code{x},
#' though the number of rows will depend on whether any subsetting (if \code{subset.row} is supplied)
#' or feature selection (if \code{use.theirs=TRUE}) was performed.
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
#' out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced))
#'
#' # make scvelo use 10 rather than the default 30 neighbors to compute moments for velocity estimation:
#' out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced),
#'               scvelo.params=list(neighbors=list(n_neighbors=10L)))
#'
#' @references
#' Bergen, V., Lange, M., Peidli, S. et al. Generalizing RNA velocity to transient cell states through dynamical modeling. Nat Biotechnol 38, 1408–1414 (2020). \url{https://doi.org/10.1038/s41587-020-0591-3}
#'
#' @author Aaron Lun, Charlotte Soneson
#' @name scvelo
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle normalizeCounts
#' @importFrom BiocSingular bsparam runPCA
#' @importFrom BiocParallel SerialParam
#' @importFrom Matrix t
.scvelo <- function(x, subset.row=NULL,
    sf.X=NULL, sf.spliced=NULL, sf.unspliced=NULL,
    use.theirs=FALSE,
    mode=c('steady_state', 'deterministic', 'stochastic', 'dynamical'),
    scvelo.params=list(), dimred=NULL, ncomponents=30, BPPARAM=SerialParam(),
    BSPARAM=bsparam())
{
    spliced <- x$spliced
    unspliced <- x$unspliced
    X <- x$X

    refdim <- as.integer(dim(spliced))
    if (!identical(refdim, as.integer(dim(unspliced))) || !identical(refdim, as.integer(dim(X)))) {
        stop("matrices in 'x' must have the same dimensions")
    }

    if (!use.theirs) {
        X <- normalizeCounts(X, sf.X, log=TRUE)
        spliced <- normalizeCounts(spliced, sf.spliced, log=FALSE)
        unspliced <- normalizeCounts(unspliced, sf.unspliced, log=FALSE)

        if (!is.null(subset.row)) {
            X <- X[subset.row,,drop=FALSE]
            spliced <- spliced[subset.row,,drop=FALSE]
            unspliced <- unspliced[subset.row,,drop=FALSE]
        }

        if (is.null(dimred)) {
            dimred <- runPCA(t(X), rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM)$x
        }
    }

    mode <- match.arg(mode)
    output <- basiliskRun(env=velo.env, fun=.run_scvelo,
        X=X, spliced=spliced, unspliced=unspliced,
        use.theirs=use.theirs, mode=mode,
        scvelo.params=scvelo.params,
        dimred=dimred, testload = c("scvelo", "anndata"))

    output
}

.run_scvelo <- function(X, spliced, unspliced, use.theirs=FALSE, mode='dynamical', scvelo.params=list(), dimred=NULL) {
    X <- t(velociraptor:::.make_np_friendly(X))
    spliced <- t(velociraptor:::.make_np_friendly(spliced))
    unspliced <- t(velociraptor:::.make_np_friendly(unspliced))

    and <- reticulate::import("anndata")
    scv <- reticulate::import("scvelo")
    sc <- reticulate::import("scanpy")
    adata <- and$AnnData(X, layers=list(spliced=spliced, unspliced=unspliced))
    adata$obs_names <- rownames(spliced)
    adata$var_names <- colnames(spliced)

    ## A supplied dimred will be used even if use.theirs=TRUE
    if (!is.null(dimred)) {
        dimred <- velociraptor:::.make_np_friendly(dimred)
        adata$obsm <- list(X_pca = dimred)
    }

    if (use.theirs) {
        do.call(scv$pp$filter_and_normalize, c(list(data=adata), scvelo.params$filter_and_normalize))
    }

    do.call(sc$pp$neighbors, c(list(adata), scvelo.params$neighbors))
    
    if (!is.null(scvelo.params$moments)) {
        if (!is.null(scvelo.params$moments$n_neighbors)) {
            stop("scvelo.params$moments$n_neighbors is deprecated; use scvelo.params$neighbors$n_neighbors instead")
        }
        if (!is.null(scvelo.params$moments$n_pcs)) {
            stop("scvelo.params$moments$n_pcs is deprecated; use scvelo.params$neighbors$n_pcs instead")
        }
    } else {
        # if unspecified, set to NULL (= None)
        # see https://github.com/theislab/scvelo/issues/1212
        scvelo.params$moments <- list(
            n_neighbors = NULL,
            n_pcs = NULL
        )
    }
    
    do.call(scv$pp$moments, c(list(data=adata), scvelo.params$moments))

    if (mode=="dynamical") {
        do.call(scv$tl$recover_dynamics, c(list(data=adata), scvelo.params$recover_dynamics))
    }

    scvelo.params$velocity$mode <- mode
    do.call(scv$tl$velocity, c(list(data=adata), scvelo.params$velocity))

    do.call(scv$tl$velocity_graph, c(list(data=adata), scvelo.params$velocity_graph))

    do.call(scv$tl$velocity_pseudotime, c(list(adata=adata), scvelo.params$velocity_pseudotime))

    if (mode=="dynamical") {
        do.call(scv$tl$latent_time, c(list(data=adata), scvelo.params$latent_time))
    }

    do.call(scv$tl$velocity_confidence, c(list(data=adata), scvelo.params$velocity_confidence))

    zellkonverter::AnnData2SCE(adata)
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
    assay.X="counts", assay.spliced="spliced", assay.unspliced="unspliced")
{
    .scvelo(list(X=assay(x, assay.X), spliced=assay(x, assay.spliced), unspliced=assay(x, assay.unspliced)), ...)
})

#' @export
#' @rdname scvelo
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
setMethod("scvelo", "SingleCellExperiment", function(x, ..., sf.X=sizeFactors(x), dimred=NULL, use.dimred=NULL) {
    if (is.null(dimred)) {
        if (!is.null(use.dimred)) {
            dimred <- reducedDim(x, use.dimred)
        }
    }
    callNextMethod(x, ..., sf.X=sf.X, dimred=dimred)
})
