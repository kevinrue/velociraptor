# This tests for the correct function of the scvelo function.
# library(testthat); library(velociraptor); source("test-scvelo.R")

library(scuttle)
sce1 <- mockSCE()
sce2 <- mockSCE()

spliced <- counts(sce1)
unspliced <- counts(sce2)

set.seed(100001)
test_that("scvelo works as expected with vanilla settings", {
    out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced))
    expect_identical(rownames(out), rownames(spliced))
    expect_identical(colnames(out), colnames(spliced))
    expect_true("velocity" %in% assayNames(out))
    expect_true("velocity_pseudotime" %in% colnames(colData(out)))

    # Throwing in some subsetting.
    subset.row <- c(50:1, 100:200)
    out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced), subset.row=subset.row)
    expect_identical(rownames(out), rownames(spliced)[subset.row])
    expect_identical(colnames(out), colnames(spliced))
    expect_true("velocity" %in% assayNames(out))
    expect_true("velocity_pseudotime" %in% colnames(colData(out)))

    # Check that the errors break.
    expect_error(scvelo(list(X=spliced[1:10,], spliced=spliced, unspliced=unspliced)), "same dimensions") 
})

set.seed(100002)
test_that("scvelo is agnostic to X when we override the dimensionality reduction", {
    dimred <- matrix(rnorm(10000), nrow=ncol(spliced))

    set.seed(0)
    ref <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced), dimred=dimred)

    set.seed(0)
    spliced2 <- spliced
    spliced2[] <- 0
    alt <- scvelo(list(X=spliced2, spliced=spliced, unspliced=unspliced), sf.X=rep(1, ncol(spliced2)), dimred=dimred)

    # TODO: find and kill that random seed!
    expect_equal(alt$velocity_pseudotime, ref$velocity_pseudotime, tol=1e-5)
})

set.seed(100003)
test_that("scvelo works correctly with SE inputs", {
    sce <- SingleCellExperiment(list(spliced=counts(sce1), unspliced=counts(sce2)))

    set.seed(0)
    ref <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced))
    set.seed(0)
    alt <- scvelo(sce, assay.X="spliced")
    expect_equal(alt$velocity_pseudotime, ref$velocity_pseudotime, tol=1e-5)

    # Behaves properly with dimreds.
    dimred <- matrix(rnorm(10000), nrow=ncol(spliced))

    set.seed(0)
    ref <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced), dimred=dimred)
    set.seed(0)
    reducedDim(sce, "PCA") <- dimred
    alt <- scvelo(sce, assay.X="spliced", use.dimred="PCA")
    expect_equal(alt$velocity_pseudotime, ref$velocity_pseudotime, tol=1e-5)
})

set.seed(100004)
test_that("scvelo works correctly with their pipeline", {
    ref <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced), use.theirs=TRUE)
    expect_s4_class(ref, "SingleCellExperiment")
})
