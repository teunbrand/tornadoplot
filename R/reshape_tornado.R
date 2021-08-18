# Reshaping ---------------------------------------------------------------

#' Reorder features of a tornado
#'
#' @param tornado A \linkS4class{TornadoExperiment} object.
#' @param assay_name A `character(1)`: one of the assay names.
#' @param bin_subset,bin_weights Define either `bin_subset` or `bin_weights` (or
#'   neither), but not both. `bin_subset` should be a `logical` or `numeric`
#'   vector with subsetting indices for the number of bins. `bin_weights` should
#'   be a `numeric` vector of weights, the same length as there are bins. If
#'   both `bin_subset` and `bin_weights` are `NULL` (default), the weights are
#'   constructed with the density function of the Laplace distribution wherein
#'   the location is centred at the middle of the bins and the scale parameter
#'   is 10% of the number of bins.
#' @param sample_subset A `character`, `logical` or `integer` vector for
#'   subsetting the samples.
#' @param decreasing A `logical(1)`: sort from high to low (`TRUE`) or from low
#'   to high (`FALSE`)?
#'
#' @return The `tornado` dataset with the features reordered.
#' @export
#' @importFrom SummarizedExperiment assayNames
#'
#' @examples
#' NULL
sort_tornado <- function(
  tornado,
  assay_name = assayNames(tornado),
  bin_subset  = NULL,
  bin_weights = NULL,
  sample_subset = NULL,
  decreasing = TRUE
) {
  assay <- match.arg(assay_name, assayNames(tornado))
  dim <- dim(tornado)

  # Resolve weights
  bin_weights <- resolve_bin_weights(dim[3], bin_subset, bin_weights)

  # Resolve sample subset
  if (length(sample_subset) < 1) {
    sample_subset <- seq_len(dim[2])
  }
  sample_subset <- as.integer(NSBS(
    sample_subset, setNames(seq_len(dim[2]), colnames(tornado))
  ))

  # Determine order
  order <- assay(tornado, assay)[, sample_subset, ]
  order <- sweep(assay(tornado, assay), 3, bin_weights, FUN = "*")
  order <- rowSums(order)
  if (decreasing) {
    order <- -order
  }

  order <- order(get_keydata(tornado, "row"), order)
  tornado[order,]
}

#' Reshape tornado to long format
#'
#' @inheritParams sort_tornado
#'
#' @return A `data.frame`
#' @export
#'
#' @examples
#' NULL
melt_tornado <- function(tornado, assay_name = assayNames(tornado)) {
  assay <- match.arg(assay_name, assayNames(tornado))

  position <- get_bin_position(tornado)
  feat_num <- get_feature_num(tornado)
  samp_dat <- get_sample_data(tornado)

  cbind.data.frame(
    samp_dat[flat_idx(tornado, 2), ],
    feat_num[flat_idx(tornado, 1), ],
    position = position[flat_idx(tornado, 3)],
    value = as.vector(assay(tornado, assay_name)),
  )
}

#' @importFrom SummarizedExperiment assay<-
norm_samples <- function(x, scale,
                         assay = assayNames(x), assay_to = assay) {
  if (length(scale) == 1) {
    scale <- rep(scale, ncol(x))
  }
  if (!is.null(names(scale)) && !is.null(colnames(x))) {
    i <- match(names(scale), colnames(x))
    scale <- scale[i]
  }
  stopifnot(
    "Cannot match up the 'scale' parameter with columns in 'x'." =
      length(scale) == ncol(x),
    "Scale cannot contain NAs" =
      all(!is.na(scale))
  )
  assay    <- match.arg(assay,    assayNames(x))
  assay_to <- match.arg(assay_to, assayNames(x))
  assay <- assay(x, assay, withDimnames = FALSE)

  assay <- sweep(assay, 2, scale, FUN = "/")
  assay(x, assay_to, withDimnames = FALSE) <- assay
  x

}

#' Flatten features
#'
#' Summarise a `TornadoExperiment` at every bin position for every sample by
#' calculating a statistic over all features in a set. The defaults measure
#' the mean, standard deviation and number of features.
#'
#' @param x A `TornadoExperiment` object.
#' @param assay An `integer(1)` or `character(1)` forwarded to the `i` argument
#'   of the [`assay()`][SummarizedExperiment] function.
#' @param measure A **named** `list` wherein each element is a `function`.
#'   Every function is expected to take a `nrow(x)` by `nbin(x)` matrix and
#'   return one of these results:
#'   \itemize{
#'    \item An atomic vector of length 1. An example of this is `nrow()`.
#'    \item An atomic vector of length `nbin(x)`. An example of this is
#'      `colMeans()`.
#'    \item An atomic matrix of dimension `nbin(x)` by `m` *with* columnames.
#'      An example of this is the [`colQuantiles`][matrixStats::colQuantiles]
#'      function.
#'   }
#'   As one might infer from the mentioned example functions, the
#'   `col*`-patterned functions are generally appropriate. Alternatively, the
#'   list elements can be `formula` with right-hand side. This follows the
#'   \pkg{rlang} [lambda syntax][rlang::as_function()].
#'
#' @return A long format `data.frame` with columns named after the `measure`
#'   argument, as well as `position`, `feature_set` and elements from
#'   `colData(x)`.
#' @export
#' @importFrom MatrixGenerics colMeans2 colSds
#'
#' @examples
#' NULL
flatten_features <- function(
  x, assay = assayNames(x),
  measure = list(mean = colMeans2, sd = colSds, n = nrow)
) {
  nbin <- nbin(x)
  measure <- resolve_measure(measure)
  assay <- match.arg(assay, assayNames(x))
  diced <- dice_tornado(x, assay = assay)

  # Collect metadata
  position <- get_bin_position(x)
  sample   <- get_sample_data(x)
  featset  <- unique(get_keydata(x, "row"))

  # Do measurements
  ans <- mapply(function(mat, i) {
    res <- lapply(measure, function(f) {
      f(mat)
    })
    res$`.index` <- i
    res$`.bin` <- seq_len(nbin)
    res
  }, mat = diced, i = seq_along(diced), SIMPLIFY = FALSE)
  ans <- resolve_listmatrix(ans, nbin)

  # Attach metadata
  ans$position <- position[ans$.bin]
  ans$feature_set <- featset[flat_idx(diced, 1)[ans$.index]]
  ans[names(sample)] <- as.list(sample[flat_idx(diced, 2)[ans$.index], ])

  # Cleanup
  ans$.bin <- ans$.index <- NULL
  ans
}

#' Collapse tornado features
#'
#' Use a summarising function to collapse a tornado along the features.
#'
#' @param tornado A `tornado_array` dataset.
#' @param fun A function to apply to a matrix build from every section of the
#'   array.
#' @param ... Extra arguments to pass to the `fun` function.
#'
#' @return A `data.frame` with `value`, `feature_set`, `position` and `sample`
#'   columns.
#' @export
#'
#' @examples
#' NULL
collapse_tornado <- function(tornado, fun = colMeans, ...) {
  tor <- tornado$tornado
  dim <- dim(tor)

  bins <- tornado$row_dat
  bins <- runValue(bins) * runLength(bins) - length(bins) * 0.5 -
    0.5 * runLength(bins)[1]

  # Split tornado along samples
  tor <- split(tor, slice.index(tor, 3))

  # Split samples along feature sets
  tor <- lapply(tor, function(x) {
    dim(x) <- dim[1:2]
    x <- split.data.frame(t(x), decode(tornado$col_dat))
    x <- lapply(x, fun)
    lens <- lengths(x)
    x <- unlist(x, use.names = FALSE)
    names <- names(lens)
    list2DF(list(
      value = x,
      feature_set = rep.int(names(lens), lens),
      position = rep(bins, length(lens))
    ))
  })
  lens <- vapply(tor, nrow, integer(1))
  tor  <- do.call(rbind, c(tor, list(make.row.names = FALSE,
                                     deparse.level = 0)))
  tor$sample <- tornado$slice_dat[rep.int(seq_along(lens), lens)]

  return(tor)
}

# Dicing ------------------------------------------------------------------

.dice_tornado <- function(x, f, ...) {
  dim <- dim(x)
  stopifnot(
    "Cannot dice tornado: 'f' must be of length 'nrow(x)'." =
      length(f) == dim[1],
    "Cannot dice array with dimensions other than 3." =
      length(dim) == 3L
  )
  f <- split(seq_along(f), f)

  i <- rep.int(seq_along(f), dim[2])
  j <- rep(seq_len(dim[2]), each = length(f))

  ans <- mapply(
    function(i, j) x[i, j, , drop = TRUE],
    i = f[i], j = j, SIMPLIFY = FALSE
  )
  dim(ans) <- c(length(f), dim[2])
  ans
}

setGeneric(
  "dice_tornado",
  function(x, f, ...) standardGeneric("dice_tornado")
)

setMethod(
  "dice_tornado", "TornadoExperiment",
  function(x, f = NULL, assay = assayNames(x)) {
    assay <- match.arg(assay, assayNames(x))
    if (missing(f) || is.null(f)) {
      f <- get_keydata(x, keytype = "row")
    }
    ans <- callGeneric(
      x = assay(x, assay),
      f = f
    )
    dimnames(ans) <- list(unique(f), get_keydata(x, keytype = "col"))
    ans
  }
)

setMethod(
  "dice_tornado", "array",
  .dice_tornado
)

# Resolve arguments -------------------------------------------------------

resolve_bin_weights <- function(nbins, bin_subset, bin_weights) {
  x <- seq_len(nbins)
  if (!is.null(bin_subset) && !is.null(bin_weights)) {
    warning("Cannot use `bin_subset` and `bin_weights` at the same time.")
    bin_subset <- NULL
  }
  # Translate bin_subset to 0/1 weights
  if (!is.null(bin_subset) && is.null(bin_weights)) {
    bin_subset  <- as.integer(NSBS(bin_subset, x))
    bin_weights <- as.numeric(x %in% bin_subset)
  }
  # Default to Laplacian density weights
  if (is.null(bin_weights)) {
    bin_range   <- range(x)
    bin_weights <- dlaplace(x, mu = 0.5 * diff(bin_range), b = 0.1 * nbins)
  }
  # Expand weights if they are of length 1
  if (length(bin_weights) == 1) {
    bin_weights <- rep(bin_weights, nbins)
  }
  # Check if length matches
  if (!is.null(bin_weights)) {
    stopifnot(
      "Bin weights do not match number of bins" =
        length(bin_weights) == nbins
    )
  }
  bin_weights
}

#' @importFrom rlang as_function is_formula
resolve_measure <- function(x, argname = "measure") {
  if (is.function(x) || is_formula(x)) {
    x <- list("func" = x)
    msg <- paste0(
      "The '", argname, "' argument was not a list but a single function. ",
      "It will have the column name 'func' in the output."
    )
    message(msg)
  }

  # Check names
  names <- names(x)

  msg <- NULL
  if (any(!nzchar(names))) {
    msg <- c(msg, paste0(
      "All functions in the '", argname, "' should be named."
    ))
  }
  if (any(c(".index", ".bin") %in% names(x))) {
    msg <- c(msg, paste0(
      "'.index' and '.bin' are prohibited names in the '",
      argname, "' argument."
    ))
  }
  # Covert formulas to functions
  formulas <- vapply(x, is_formula, logical(1))
  if (any(formulas)) {
    x[formulas] <- lapply(x[formulas], as_function)
  }
  # Check all are functions
  is_fun <- vapply(x, is.function, logical(1))
  if (any(!is_fun)) {
    msg <- c(msg, paste0(
      "The '", argname, "' argument should be a list of functions."
    ))
  }

  if (length(msg)) {
    stop(paste0(msg, collapse = "\n"), call. = FALSE)
  }

  return(x)
}

resolve_listmatrix <- function(x, nbin) {
  nans <- unique(lengths(x))
  stopifnot(length(nans) == 1)
  dim <- c(nans, length(x))
  nms <- names(x[[1]])

  # Bind columns
  if (dim[2] > 1) {
    x <- with_dim(unlist(x, FALSE, FALSE), dim)
    x <- apply(x, 1, function(x) {bindROWS(x[[1]], x[-1])})
  } else {
    x <- unlist(x, FALSE, FALSE)
  }
  x <- setNames(x, nms)

  # Resolve rows
  rows <- vapply(x, NROW, integer(1))
  x[rows == dim[2]] <- lapply(x[rows == dim[2]], rep, each = nbin)
  delete <- vapply(x, NROW, integer(1)) != (dim[2] * nbin)
  if (any(delete)) {
    msg <- names(delete)[delete]
    if (length(msg) > 1) {
      msg <- paste0(paste0(msg[-length(msg)], collapse = ", "),
                    " and ", msg[length(msg)])
    }
    msg <- paste0("The following measurements had inappropriate lengths: ",
                  msg)
    warning(msg, call. = FALSE)
    x[delete] <- NULL
  }

  # Resolve columns
  cols <- vapply(x, NCOL, integer(1))
  alter <- names(x)[cols != 1L]
  alter <- lapply(alter, function(i) {
    out <- as.list(as.data.frame(x[[i]]))
    outnames <- names(out) %||% seq_len(ncol(out))
    setNames(out, paste0(i, "_", outnames))
  })
  alter <- unlist(alter, FALSE)
  x <- c(x[cols == 1L], alter)
  list2DF(x)
}

# Helpers -----------------------------------------------------------------

# Based on https://en.wikipedia.org/wiki/Laplace_distribution
dlaplace <- function(x, mu = 0, b = 1) {
  (1 / (2 * b)) * exp(-1 * (abs(x - mu)/b))
}

get_bin_position <- function(x, extremes = FALSE) {
  pos <- binData(x)
  if ("range" %in% colnames(pos) && inherits(pos$range, "Ranges")) {
    pos <- (start(pos$range) + end(pos$range)) / 2
  } else {
    pos <- as.vector(pos[[binKey(x)]])
    if (!is.numeric(pos)) {
      pos <- match(pos, sort(unique(pos)))
    }
  }
  pos
}

get_feature_num <- function(x) {
  feat <- get_keydata(x, "row")
  nr <- unlist(lapply(split(feat, feat), seq_along), FALSE, FALSE)
  data.frame(
    feature_set = as.vector(feat),
    feature_nr  = as.vector(nr)
  )
}

get_sample_data <- function(x) {
  # Put key as first column
  key <- colKey(x)
  x <- colData(x)
  notkey <- setdiff(colnames(x), key)
  x <- x[, c(key, notkey)]
  as.data.frame(x)
}

with_dim <- function(x, dim) {
  dim(x) <- dim
  x
}

bind_list <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  bindROWS(x[[1]], x[-1])
}

flat_idx <- function(x, margin = 1) {
  as.vector(slice.index(x, margin))
}
