#' @importFrom DelayedArray is_sparse
.make_np_friendly <- function(x) {
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
    } else {
        as.matrix(x) 
    }
}

.extractor_python_dict <- function(thing, names, single=FALSE) {
    if (single) {
        values <- lapply(names, function(x) thing[x])
    } else {
        values <- lapply(names, function(x) thing[[x]])
    }
    names(values) <- names
    values
}
