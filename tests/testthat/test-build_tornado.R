test_that("tornados can be build from granges", {
  f <- dummy_features()
  d <- dummy_granges_data()

  tor <- build_tornado(f, d, width = 3000)

  expect_equal(dim(tor), c(2, 1, 120))

  d <- GenomicRanges::GRangesList("A" = d, "B" = d)

  tor <- build_tornado(f, d, width = 3000)

  expect_equal(dim(tor), c(2, 2, 120))
})

test_that("tornados can be build from BigWigFiles", {
  f <- dummy_features()
  d <- dummy_bigwig()

  tor <- build_tornado(f, d, width = 3000)

  expect_equal(dim(tor), c(2, 1, 120))

  d <- rtracklayer::BigWigFileList(
    setNames(rep(path(d), 2), c("A", "B"))
  )

  tor <- build_tornado(f, d, width = 3000)

  expect_equal(dim(tor), c(2, 2, 120))

  d <- path(d)

  tor <- build_tornado(f, d, width = 3000)

  expect_equal(dim(tor), c(2, 2, 120))

  unlink(d)
})

test_that("tornados can be build from TabixFiles", {

  f <- dummy_features()
  d <- dummy_tabix()

  tor <- build_tornado(f, d, width = 3000)

  expect_equal(dim(tor), c(2, 1, 120))

  tor <- build_tornado(f, d, width = 3000,
                       barcode_groups = list("A" = "A", "B" = "B"),
                       barcode_column = 4)

  expect_equal(dim(tor), c(2, 2, 120))

})

test_that("tornados flip reverse strand", {

  f <- dummy_features()
  strand(f) <- c("+", "-")
  f <- resize(f, width(f) - 200)
  d <- dummy_granges_data()

  tor <- build_tornado(f, d, width = 2000)
  ctl <- build_tornado(unstrand(f), d, width = 2000)
  tor <- assay(tor)
  ctl <- assay(ctl)
  expect_equal(tor[1,,], ctl[1,,])
  expect_equal(tor[2,,], rev(ctl[2,,]))
  expect_false(all(tor[2,,] == ctl[2,,]))
})

test_that("tornados can handle weird seqlenghts", {

  f <- dummy_features()
  f <- GenomicRanges::shift(f, 500)

  d <- dummy_granges_data()

  tor <- build_tornado(f, d, width = 3000)
  tor <- assay(tor)
  expect_equal(sum(colSums(tor)[101:120]), 0)
})

test_that("check_seqlevels  throws appropriate errors", {

  expect_error(
    check_seqlevels(c("chr1", "chr2"), c("1", "2")),
    "The features and data have no common sequences."
  )

  expect_error(
    check_seqlevels(c("chr1", "chr2"), c("chr3", "chr4")),
    "The features and data have no common sequences."
  )

  expect_warning(
    check_seqlevels(c("chr1", "chr2"), c("chr1", "chr3")),
    "The following sequence names were missing"
  )

})

test_that("bins can be resolved", {
  x <- resolve_bins(2000, binwidth = 25)
  expect_equal(x, Rle(1:80, 25))

  expect_error(resolve_bins(6, binwidth = 4),
               "`width` is not a multiple of binwidth")

  expect_error(resolve_bins(5),
               "Define 2 out of 3")

  x <- resolve_bins(12, 3, 3) # binwidth gets overruled
  expect_equal(x, Rle(1:3, 4))

  x <- resolve_bins(binwidth = 4, nbin = 3)
  expect_equal(x, Rle(1:3, 4))

  expect_error(resolve_bins(12, nbin = 5),
               "`width` does not fit an integer of `nbin`")

})

test_that("features can be resolved", {

  bins <- resolve_bins(2000, 25)

  expect_error(
    resolve_features("A", bins),
    "be a `GRanges` or `GRangesList` object"
  )

  f <- rep(GRanges("chr1:3000-3100"), 50001)
  expect_warning(
    resolve_features(f, bins),
    "consider downsampling the features"
  )

  f <- GRanges(c("chr1:100-200", "chr1:5000-6000"))
  expect_warning(
    resolve_features(f, bins),
    "due to negative start sizes"
  )


})
