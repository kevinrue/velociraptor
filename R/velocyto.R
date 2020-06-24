#' Velocity calculations with \pkg{velocyto}
#'
#' Performs RNA velocity calculations using the original implementation
#' of the steady-state model by La Manno et al. (2018).
#'
#' @param x A list of count matrices of the same dimensions, with genes in rows and cells in samples.
#' This should contain \code{"spliced"} and \code{"unspliced"} (and possibly \code{"ambiguous"}, see Details).
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing three such matrices among its assays.
#' @inheritParams scvelo
#' @param use.theirs Logical scalar indicating whether \pkg{velocyto}'s gene filtering and normalization should be used.
#' @param velocyto.params List of lists containing arguments for individual \pkg{velocyto} functions, see details below.
#' @param assay.ambiguous  An integer scalar or string specifying the assay \code{x} containing the ambiguous counts.
#' If \code{NULL}, no attempt is made to extract this count matrix, see Details.
#' 
#' @details
#' This function uses the \pkg{velocyto} Python package (\url{https://pypi.org/project/velocyto/}) for RNA velocity calculations.
#' All matrices are saved to a temporary \code{.loom} file that is subsequently read by \pkg{velocyto} functions.
#' This is necessary as \pkg{velocyto} itself does not provide an entry point for matrices that are already in memory.
#'
#' For consistency with other Bioconductor workflows, we perform as many standard steps in R as we can
#' before starting the velocity calculations with \pkg{scVelo}.
#' This involves:
#' \enumerate{
#' \item Size factor calculations, using \code{\link{librarySizeFactors}} on entries of \code{x}.
#' If \code{sf.*} arguments are supplied, these are used directly.
#' \item Subsetting all matrices to \code{subset.row}, most typically to a subset of interest, e.g., highly variable genes.
#' Note that, if set, any subsetting is done \emph{after} size factor calculations so that library sizes are correctly computed.
#' In addition, the \code{.loom} file will only contain the specified subset to reduce the I/O overhead.
#' }
#'
#' Users can set \code{use.theirs=TRUE} to directly use the entire \pkg{velocyto} normalization and filtering pipeline.
#' This ignores all of the size factor arguments (\code{sf.*}) as well as \code{subset.row}.
#'
#' All individual \pkg{velocyto} function calls are performed with their default values,
#' but this can be changed by including named parameter lists in \code{velocyto.params}.
#' Supported arguments are:
#' \itemize{
#' \item \code{score_cv_vs_mean}, for CV-mean calculations if \code{use.theirs=TRUE}.
#' \item \code{filter_genes}, for feature selection if \code{use.theirs=TRUE}.
#' \item \code{normalize}, for normalization if \code{use.theirs=TRUE}.
#' \item \code{perform_PCA}
#' \item \code{knn_imputation}
#' \item \code{fit_gammas}
#' \item \code{predict_U}
#' \item \code{calculate_velocity}
#' \item \code{calculate_shift}
#' \item \code{extrapolate_cell_at_t}
#' }
#'
#' It does not seem that the \code{"ambiguous"} entry in \code{x} is used for anything,
#' but it is nonetheless necessary to supply \emph{some} value as \pkg{velocyto} will fail otherwise.
#' If such an entry is not present in \code{x}, an empty sparse matrix will be created instead.
#'
#' Upon first use, this function will instantiate a conda environment containing the \pkg{velocyto} package.
#' This is done via the \pkg{basilisk} package - see the documentation for that package for trouble-shooting.
#'
#' @return
#' A \linkS4class{SummarizedExperiment} is returned containing the output of the velocity calculations.
#' TODO: SOME DOCUMENTATION REQUIRED.
#'
#' @references
#' La Manno G, Soldatov R, Zeisel A et al. (2018).
#' RNA velocity of single cells.
#' \emph{Nature} 560, 494-498
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
#' out <- velocyto(list(spliced=spliced, unspliced=unspliced))
#'
#' @name velocyto
NULL

#' @importFrom scuttle librarySizeFactors
#' @importFrom Matrix sparseMatrix
.velocyto <- function(x, subset.row=NULL,
    sf.spliced=NULL, sf.unspliced=NULL,
    use.theirs=FALSE,
    velocyto.params=list())
{
    if (!use.theirs) {
        if (is.null(sf.spliced)) {
            sf.spliced <- librarySizeFactors(x$spliced)
        }
        if (is.null(sf.unspliced)) {
            sf.unspliced <- librarySizeFactors(x$unspliced)
        }

        if (!is.null(subset.row)) {
            x <- lapply(x, function(x) x[subset.row,,drop=FALSE])
        }
    }

    if (!"ambiguous" %in% names(x)) {
        x$ambiguous <- sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=dim(x$spliced))
    }

    x0 <- c(list(matrix=x$spliced), x) # adding an extra matrix to avoid loss of the 'spliced' layer.
    scle <- LoomExperiment::SingleCellLoomExperiment(x0)
    path <- tempfile(fileext=".loom")
    on.exit(unlink(path))
    LoomExperiment::export(scle, con=path)

    output <- basiliskRun(env=velocyto.env, fun=.run_velocyto,
        path=path, use.theirs=use.theirs, velocyto.params=velocyto.params,
        sf.spliced=sf.spliced, sf.unspliced=sf.unspliced)
}

#' @importFrom reticulate import
.run_velocyto <- function(path, use.theirs=FALSE, velocyto.params=list(), sf.spliced, sf.unspliced) {
    vcy <- import("velocyto")
    vdata <- vcy$VelocytoLoom(path)

    if (use.theirs) {
        do.call(vdata$score_cv_vs_mean, c(list(), velocyto.params$score_cv_vs_mean))
        do.call(vdata$filter_genes, c(list(), velocyto.params$filter_genes))
    }

    if (use.theirs) {
        do.call(vdata$normalize, velocyto.params$normalize)
    } else {
        vdata$normalize(which="S", relative_size=sf.spliced)
        vdata$normalize(which="U", relative_size=sf.unspliced)
    }

    do.call(vdata$perform_PCA, c(list(), velocyto.params$perform_PCA))
    do.call(vdata$knn_imputation, c(list(), velocyto.params$knn_imputation))

    do.call(vdata$fit_gammas, c(list(), velocyto.params$fit_gammas))
    do.call(vdata$predict_U, c(list(), velocyto.params$predict_U))
    do.call(vdata$calculate_velocity, c(list(), velocyto.params$calculate_velocity))
    do.call(vdata$calculate_shift, c(list(), velocyto.params$calculate_shift))
    do.call(vdata$extrapolate_cell_at_t, c(list(), velocyto.params$extrapolate_cell_at_t))

    SummarizedExperiment(list(velocity=vdata$velocity, extrapolated=vdata$Sx_sz_t))
}

#' @export
#' @rdname velocyto
setGeneric("velocyto", function(x, ...) standardGeneric("velocyto"))

#' @export
#' @rdname velocyto
setMethod("velocyto", "ANY", .velocyto)

#' @export
#' @rdname velocyto
setMethod("velocyto", "SummarizedExperiment", function(x, ...,
    assay.spliced="spliced", assay.unspliced="unspliced", assay.ambiguous=NULL)
{
    stuff <- list(spliced=assay(x, assay.spliced), unspliced=assay(x, assay.unspliced))
    if (is.null(assay.ambiguous)) {
        stuff$ambiguous <- assay(x, assay.ambiguous)
    }
    .velocyto(stuff, ...)
})
