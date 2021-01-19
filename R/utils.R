#' @importFrom DelayedArray is_sparse
.make_np_friendly <- function(x) {
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
    } else {
        as.matrix(x)
    }
}

#' @importFrom DelayedArray t
.extractor_python_dict <- function(thing, names, single=FALSE, transpose=FALSE) {
    if (single) {
        values <- lapply(names, function(x) {
            if (transpose) t(thing[x])
            else thing[x]
        })
    } else {
        values <- lapply(names, function(x) {
            if (transpose) t(thing[[x]])
            else thing[[x]]
        })
    }
    names(values) <- names
    values
}
