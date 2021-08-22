# Reshaping ---------------------------------------------------------------

#' Reorder features
#'
#' The `sort_tornado()` function reorders the features of a
#' \linkS4class{TornadoExperiment} object based on the weighted sum across all
#' bins.
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
#' @param ... Not currently in use.
#'
#' @return The `tornado` dataset with the features reordered.
#' @export
#' @importFrom SummarizedExperiment assayNames
#' @family tornado utilities
#'
#' @examples
#' # A tornado that isn't sorted
#' tor <- dummy_tornado()[c(2,1,4,3),]
#'
#' # Sorting based on a particular sample
#' x <- sort_tornado(tor, sample_subset = "dummy_ctrl")
#'
#' # Sorting based on a subset of bins
#' x <- sort_tornado(tor, bin_subset = 21:30)
#'
#' # Sorting with custom weights
#' nbin <- nbin(tor)
#' w <- dnorm(seq_len(nbin), mean = (nbin + 1)/2, sd = 5)
#' x <- sort_tornado(tor, bin_weights = w)
sort_tornado <- function(
  tornado,
  assay_name = assayNames(tornado),
  bin_subset  = NULL,
  bin_weights = NULL,
  sample_subset = NULL,
  decreasing = TRUE,
  ...
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
  metadata(tornado) <- c(metadata(tornado), list(sorted = TRUE))
  tornado[order,]
}

#' Reshape tornado to long format
#'
#' The `melt_tornado` emulates for the \linkS4class{TornadoExperiment} class
#' what the `melt()` function does in the \pkg{reshape2} package. It converts
#' the data to long format, which is the format of choice in the \pkg{tidyverse}
#' family of packages.
#'
#' @inheritParams sort_tornado
#'
#' @return A `data.frame` with `prod(dim(tornado))` rows.
#'
#' @note For plotting purposes, the [`prep_tornado()`][prep_tornado()] function
#' does a more efficient job than `melt_tornado()`.
#'
#' @export
#' @family tornado utilities
#'
#' @examples
#' # Melt a tornado
#' tor <- dummy_tornado()
#' # Reshaping to long format
#' df  <- melt_tornado(tor)
#'
#' # Building a tornado plot from scratch
#' require(ggplot2)
#' ggplot(df, aes(position, feature_nr)) +
#'   geom_raster(aes(fill = value)) +
#'   facet_grid(feature_set ~ sample_name)
melt_tornado <- function(tornado, assay_name = assayNames(tornado)) {
  assay <- match.arg(assay_name, assayNames(tornado))

  position <- get_bin_position(tornado)
  feat_num <- get_feature_num(tornado)
  samp_dat <- get_sample_data(tornado)

  ans <- cbind.data.frame(
    samp_dat[flat_idx(tornado, 2), ],
    feat_num[flat_idx(tornado, 1), ],
    position = position[flat_idx(tornado, 3)],
    value = as.vector(assay(tornado, assay_name))
  )
  rownames(ans) <- NULL
  ans
}


#' Normalise samples
#'
#' Normalise the samples in a \linkS4class{TornadoExperiment} by dividing their
#' values by some scaling factor.
#'
#' @inheritParams sort_tornado
#' @param scale A `numeric` vector or length `ncol(tornado)` with normalisation
#'   factors.
#' @param assay_to A `character(1)` for the target assay name for the normalised
#'   data. Defaults to the `assay_name` assay, thereby replacing input data.
#'
#' @return A \linkS4class{TornadoExperiment} object.
#' @export
#' @family tornado utilities
#' @importFrom SummarizedExperiment assay<-
#'
#' @examples
#' tor <- dummy_tornado()
#' tor <- norm_tornado(tor, c(1, 2))
norm_tornado <- function(
  tornado, scale,
  assay_name = assayNames(tornado), assay_to = assay_name) {
  if (length(scale) == 1) {
    scale <- rep(scale, ncol(tornado))
  }
  if (!is.null(names(scale)) && !is.null(colnames(tornado))) {
    i <- match(names(scale), colnames(tornado))
    scale <- scale[i]
  }
  stopifnot(
    "Cannot match up the 'scale' parameter with columns in 'x'." =
      length(scale) == ncol(tornado),
    "Scale cannot contain NAs" =
      all(!is.na(scale))
  )
  assay_name <- match.arg(assay_name, assayNames(tornado))
  assay <- assay(tornado, assay_name, withDimnames = FALSE)

  assay <- sweep(assay, 2, scale, FUN = "/")
  assay(tornado, assay_to, withDimnames = FALSE) <- assay
  tornado

}

#' Flatten features
#'
#' Summarise a `TornadoExperiment` at every bin position for every sample by
#' calculating a statistic over all features in a set. The defaults measure
#' the mean, standard deviation and number of features.
#'
#' @inheritParams sort_tornado
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
#' @family tornado utilities
#' @importFrom MatrixGenerics colMeans2 colSds
#'
#' @examples
#' # Summarise features
#' tor <- dummy_tornado()
#' df  <- flatten_features(tor)
#'
#' # Estimate standard error of the mean
#' df$se <- sqrt(df$sd^2 / df$n)
#'
#' # Plotting the result
#' require(ggplot2)
#' ggplot(df, aes(position, fill = feature_set)) +
#'   geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.3) +
#'   geom_line(aes(y = mean, colour = feature_set)) +
#'   facet_wrap(~ sample_name)
#'
#' # Calculating alternative metrics
#' require(matrixStats)
#' measure <- list(median = matrixStats::colMedians,
#'                 mad = matrixStats::colMads,
#'                 n = nrow)
#' df <- flatten_features(tor, measure = measure)
flatten_features <- function(
  tornado, assay_name = assayNames(tornado),
  measure = list(mean = colMeans2, sd = colSds, n = nrow)
) {
  nbin <- nbin(tornado)
  measure <- resolve_measure(measure)
  assay_name <- match.arg(assay_name, assayNames(tornado))
  diced <- dice_tornado(tornado, assay = assay_name)

  # Collect metadata
  position <- get_bin_position(tornado)
  sample   <- get_sample_data(tornado)
  featset  <- unique(get_keydata(tornado, "row"))

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
    msg <- comma_and(names(delete)[delete])
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

flat_idx <- function(x, margin = 1) {
  as.vector(slice.index(x, margin))
}
