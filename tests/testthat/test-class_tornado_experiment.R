

test_that("Constructors make valid objects", {
  x <- new_TornadoExperiment()
  expect_true(validObject(x))
  x <- TornadoExperiment()
  expect_true(validObject(x))

  x <- substitute(new_TornadoExperiment(
    binKey = 1
  ))
  expect_error(eval(x), 'got class "numeric", should be or extend class "character"')

  x <- substitute(TornadoExperiment(
    assays = SimpleList(w = array(1:24, 2:4)),
    rowRanges = GRanges(c("chr1:100-200", "chr2:100-200")),
    colData = DataFrame(y = 1:3),
    binData = DataFrame(z = 1:10)
  ))
  expect_error(eval(x), "nb of bins")

  tor <- dummy_tornado()
  expect_output(print(tor), "binData names\\(2\\): bin_id range")
})

test_that("Getters work", {

  tor <- dummy_tornado()

  expect_equal(dim(tor), c(4, 2, 50))

  expect_s4_class(binData(tor), "DFrame")
  expect_equal(dim(binData(tor)), c(50, 2))

  expect_s4_class(colData(tor), "DFrame")
  expect_equal(dim(colData(tor)), c(2, 2))

  expect_s4_class(rowData(tor), "DFrame")
  expect_equal(dim(rowData(tor)), c(4, 2))

  expect_null(binnames(tor))
  expect_equal(nbin(tor), 50)

  expect_equal(dimnames(tor), list(c("A", "A", "B", "B"),
                                   c("dummy_ctrl", "dummy_treat"), NULL))

  expect_equal(colKey(tor), "sample_name")
  expect_equal(rowKey(tor), "set")
  expect_equal(binKey(tor), "bin_id")
  expect_equal(get_keys(tor), setNames(c("set", "sample_name", "bin_id"),
                                       c("row", "col", "bin")))

})

test_that("Setters work", {

  # binData
  tor <- dummy_tornado()
  b <- binData(tor)
  expect_null(rownames(b))
  b$new_column <- paste0("A", seq_len(nbin(tor)))
  rownames(b) <- b$new_column
  expect_error(
    {binData(tor) <- b[1:10,]},
    "must equal nbin of object"
  )
  binData(tor) <- b
  expect_equal(binnames(tor), b$new_column)

  # dimnames
  dimnames(tor) <- NULL
  expect_null(binnames(tor))
  x <- assay(tor)
  expect_equal(dimnames(assay(tor)),
               list(NULL, NULL, NULL))
  expect_equal(dimnames(tor), NULL)

  # binnames
  binnames(tor) <- b$new_column
  expect_equal(lengths(dimnames(tor)), c(0, 0, 50))

  colKey(tor) <- "argument"
  expect_equal(tor@colKey, "argument")
  expect_error({colKey(tor) <- "nonsense"}, "does not exist")

  rowKey(tor) <- "argument"
  expect_equal(tor@rowKey, "argument")
  expect_error({rowKey(tor) <- "nonsense"}, "does not exist")

  binKey(tor) <- "new_column"
  expect_equal(tor@binKey, "new_column")
  expect_error({binKey(tor) <- "nonsense"}, "does not exist")

})


test_that("tornado can be subsetted in all three dimensions", {

  tor <- dummy_tornado()

  expect_equal(dim(tor), c(4, 2, 50))

  expect_equal(dim(tor[1:2,]), c(2, 2, 50))

  expect_equal(dim(tor[, 1]), c(4, 1, 50))

  dor <- tor[, , 11:25]
  expect_equal(dim(dor), c(4, 2, 15))

  expect_equal(dim(binData(dor)), c(15, 2))
})

