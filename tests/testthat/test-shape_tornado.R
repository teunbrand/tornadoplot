
test_that("tornado can be sorted", {
  tor <- dummy_tornado()
  alt <- tor[c(2,1,4,3),]
  sor <- sort_tornado(alt)
  expect_equal(assay(tor), assay(sor))
  expect_equal(assay(sor), assay(alt)[c(2,1,4,3),,])

  sor <- sort_tornado(tor, decreasing = FALSE)
  expect_equal(assay(sor), assay(alt))
})

test_that("bin weigths are correct", {
  nbin <- 50

  w <- resolve_bin_weights(nbin, 1:25, NULL)
  expect_equal(w, c(rep(1, 25), rep(0, 25)))

  w <- resolve_bin_weights(nbin, NULL, 1)
  expect_equal(w, rep(1, 50))

  w <- resolve_bin_weights(nbin, NULL, NULL)
  expect_equal(which.max(w), 24)

  expect_warning(
    resolve_bin_weights(nbin, 1:25, dnorm(1:50, mean = 25, sd = 10)),
    "Cannot use `bin_subset` and `bin_weights` at the same time."
  )
})

test_that("tornado can be melted", {
  tor <- dummy_tornado()

  df <- melt_tornado(tor)

  expect_equal(nrow(df), prod(dim(tor)))
  expect_equal(sum(is.na(df$value)), 0)
})

test_that("tornado can be normalised", {
  tor <- dummy_tornado()

  nor <- norm_tornado(tor, scale = c(0.5, 1))
  nor <- assay(nor)

  expect_equal(nor[,1,], nor[,2,])

  nor <- norm_tornado(tor, scale = c(dummy_treat = 1, dummy_ctrl = 0.5))
  nor <- assay(nor)
  expect_equal(nor[,1,], nor[,2,])

  nor <- norm_tornado(tor, scale = 1)
  expect_equal(assay(tor), assay(nor))
})

test_that("tornado can be flattened", {
  tor <- dummy_tornado()

  df <- flatten_features(tor)

  expect_equal(nrow(df), prod(2, ncol(tor), nbin(tor)))


  expect_warning(
    flatten_features(tor, measure = list(scale = scale)),
    "inappropriate lengths: 'scale'"
  )
})

test_that("flatten tornado resolves issues in measurement", {

  f <- resolve_measure(list(x = ~ nrow(.x)))
  expect_s3_class(f$x, "rlang_lambda_function")

  expect_message(resolve_measure(nrow), "was not a list but a single function")

  expect_error(
    resolve_measure(list(x = nrow, ncol)),
    "should be named"
  )
  expect_error(
    resolve_measure(list(.index = nrow)),
    "prohibited names"
  )
  expect_error(
    resolve_measure(list(x = "A")),
    "should be a list of functions"
  )

})

