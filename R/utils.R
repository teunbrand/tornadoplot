# Argument checking -------------------------------------------------------

#' Argument checking
#'
#' These functions are for checking that arguments of parent functions are as
#' expected.
#'
#' @param x An argument to check.
#' @param type A vector of the expected type, typically of length 0.
#' @param length An `integer(1)` or `NULL` for the expected length. When set to
#'   `NULL`, length can be anything.
#' @param allow_NAs A `logical(1)`, when `TRUE`, allows the `x` argument to have
#'   `NAs`. If `FALSE`, raises error when `x` contains `NA`s.
#' @param arg A `character(1)` for the argument being checked. This is used to
#'   identify the argument in error messaging. Defaults to the symbol used for
#'   `x`.
#' @param ... Not currently used.
#'
#' @return It returns `x` if it satisfies the criteria, or `x` coerced to the
#'   `type` class when appropriate.
#' @name norm_type
#'
#' @examples
#' # Rounding of `numeric` to `integer`
#' norm_type(1.26, integer())
norm_type_or_null <- function(x, type, length = 1, allow_NAs = FALSE,
                              arg = deparse(substitute(x))) {
  force(arg)
  if (is.null(x)) {
    return(NULL)
  } else {
    norm_type(x, type, length = length, allow_NAs = allow_NAs,
              arg = arg)
  }
}

# Type expectations -------------------------------------------------------

#' @export
#' @rdname norm_type
setGeneric(
  "norm_type", signature = c("x", "type"),
  function(x, type, length = 1, allow_NAs = FALSE,
           arg = deparse(substitute(x))) {
    standardGeneric("norm_type")
  }
)

#' @export
#' @rdname norm_type
setMethod(
  "norm_type", signature(x = "ANY", type = "numeric"),
  function(x, type, length = 1, allow_NAs = FALSE,
           arg = deparse(substitute(x))) {
    force(arg)
    x <- expect_length(x, length, arg)
    x <- expect_NAs(x, allow_NAs, arg)
    if (!is.numeric(x)) {
      msg <- paste0("`", arg, "` cannot be interpreted as ",
                    class(type)[[1]], ".")
      x <- tryCatch(
        {as.numeric(x)},
        error   = function(cond) {stop(msg, call. = FALSE); message(cond)},
        warning = function(cond) {stop(msg, call. = FALSE); message(cond)}
      )
    }
    return(x)
  }
)

#' @export
#' @rdname norm_type
setMethod(
  "norm_type", signature(x = "ANY", type = "integer"),
  function(x, type, length = 1, allow_NAs = FALSE,
           arg = deparse(substitute(x))) {
    force(arg)
    x <- norm_type(x, numeric(), length = length, allow_NAs = allow_NAs,
                   arg = arg)
    if (is.double(x)) {
      x <- as.integer(round(x, 0))
    }
    x
  }
)

# Length / validity expectations ------------------------------------------

expect_length <- function(x, length, arg = deparse(substitute(x))) {
  force(arg)
  if (is.null(length)) {
    return(x)
  }
  if (length(x) != length) {
    msg <- paste0("`", arg, "` has length ", length(x), ", expected length",
                  " is ", length, ".")
    stop(msg, call. = FALSE)
  }
  x
}

expect_NAs <- function(x, allow_NAs = FALSE, arg = deparse(substitute(x))) {
  force(arg)
  if (anyNA(x) && !allow_NAs) {
    msg <- paste0("`", arg, "` contains NAs.")
    stop(msg, call. = FALSE)
  }
  x
}

# Internals ---------------------------------------------------------------

.grab_internals <- function() {
  objects <- matrix(c(
    "replaceSlots", "BiocGenerics",
    "extract_Nindex_from_syscall", "DelayedArray",
    "normarg_names", "S4Vectors"
  ), ncol = 2, byrow = TRUE)
  funs <- setNames(nm = objects[, 1])
  out  <- mapply(function(f, ns) {
    getFromNamespace(f, ns)
  }, f = funs, ns = objects[, 2], SIMPLIFY = FALSE)
}

.int <- .grab_internals()

# Conveniences ------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

comma_and <- function(x, sep = ", ", last = " and ", quote = "'") {
  x <- paste0(quote, x, quote)
  if (length(x) < 2) {
    return(x)
  }
  paste0(
    paste0(x[-length(x)], collapse = sep),
    last, x[length(x)]
  )
}
