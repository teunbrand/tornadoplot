

test_that("tornado_plot works", {
  f <- dummy_features()
  d <- dummy_granges_data()

  grDevices::png({tmp <- tempfile(fileext = ".png")}, 100, 100)
  x <- tornado_plot(d, features = f, width = 3000)
  grDevices::dev.off()
  expect_true(file.exists(tmp))
  unlink(tmp)

  expect_s4_class(x, "TornadoExperiment")
  expect_equal(dim(x), c(2, 1, 120))

  expect_s3_class(metadata(x)$plot, "ggplot")

  cd <- colData(x)
  expect_equal(cd$sample_name, "data")
  expect_equal(cd$argument, "data")

  fd <- rowData(x)
  expect_equal(fd$set, Rle(1, 2))
  expect_equal(fd$argument, Rle("f", 2))

  bd <- binData(x)
  expect_equal(bd$bin_id, 1:120)
  expect_equal(start(bd$range), seq(-1499, 1476, by = 25))
  expect_equal(end(bd$range), seq(-1475, 1500, by = 25))

  x <- prep_tornado(x)
  expect_s3_class(x, "tornado_df")
  expect_s3_class(x$tornado, "tornado_list")

  expect_equal(as.numeric(x[, 1:4]), c(-1499.5, 1500.5, 0.5, 2.5))
})
