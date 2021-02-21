# This tests the embedding functionality.
# library(testthat); library(velociraptor); source("test-embed.R")

library(scuttle)
sce1 <- mockSCE()
sce2 <- mockSCE()

spliced <- counts(sce1)
unspliced <- counts(sce2)

set.seed(100001)
out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced))

test_that("embedVelocity works correctly", {
    tsne.results <- matrix(rnorm(2*ncol(out)), ncol=2)
    projected <- embedVelocity(tsne.results, out)
    expect_identical(dim(tsne.results), dim(projected))

    # Same results inside an SCE.
    reducedDim(sce1, "TSNE") <- tsne.results
    projected2 <- embedVelocity(sce1, out, use.dimred="TSNE")
    expect_identical(projected, projected2)
})

test_that("gridVectors works correctly", {
    tsne.results <- matrix(rnorm(2*ncol(out)), ncol=2)
    projected <- embedVelocity(tsne.results, out)

    out <- gridVectors(tsne.results, projected)
    expect_identical(ncol(out), 4L)

    expect_warning(out <- gridVectors(tsne.results, projected, return.intermediates=TRUE))
    expect_type(out, "list")
    expect_length(out, 7L)
    
    out <- gridVectors(tsne.results, projected, as.data.frame=FALSE)
    expect_type(out, "list")

    # Same results inside an SCE.
    reducedDim(sce1, "TSNE") <- tsne.results
    out2 <- gridVectors(sce1, projected, use.dimred="TSNE", as.data.frame=FALSE)
    expect_identical(out, out2)
})
