
test_that("TornadoExperiments can be autoplotted", {
  tor <- dummy_tornado()

  g <- autoplot(tor)
  expect_s3_class(g, "ggplot")

})


test_that("Autoplot/layer methods throw error if varnames not present", {

  df <- prep_tornado(dummy_tornado())

  df2 <- df
  colnames(df2)[1:2] <- c("zmin", "zmax")

  expect_error(autolayer(df2),
               "Cannot find all relevant variables.")

  df2 <- df
  colnames(df2)[6] <- "zample"

  expect_message(
    facet_tornado(df2),
    "Don't know what should be considered a sample"
  )

  f1 <- facet_tornado(df)[[1]]
  f2 <- facet_tornado(df[1:2,])[[1]]
  expect_equal(f1$params$free, list(x = FALSE, y = TRUE))
  expect_equal(f2$params$free, list(x = FALSE, y = FALSE))

})
