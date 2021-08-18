# Class definitions -------------------------------------------------------

#' @rdname TornadoExperiment
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
new_TornadoExperiment <- setClass(
  "TornadoExperiment",
  contains = "RangedSummarizedExperiment",
  representation(
    binData = "DataFrame",
    colKey  = "character",
    rowKey  = "character",
    binKey  = "character"
  ),
  prototype(
    binData = new("DFrame"),
    colKey  = character(),
    rowKey  = character(),
    binKey  = character()
  )
)

# Constructor -------------------------------------------------------------

#' TornadoExperiment objects
#'
#' The TornadoExperiment class extends the \linkS4class{SummarizedExperiment}
#' class, inheriting their methods and providing greater control over a
#' 3rd dimension of assay data. This constructor is not intended to be called
#' directly and one should see the TornadoExperiment class merely as a data
#' container class. Instead, TornadoExperiment objects should spawn from the
#' use of the [`build_tornado()`][build_tornado()] function.
#'
#' @inheritDotParams SummarizedExperiment::SummarizedExperiment
#' @param binData A \linkS4class{DataFrame} object with descriptions of the
#'   bins in the 3rd dimension of the array.
#' @param colKey,rowKey,binKey A `character(1)` containing a column name in the
#'  corresponding data that identify samples, feature sets or bins respectively.
#'  These keys exist mainly such that plotting methods can extract relevant data
#'  less ambiguously.
#'
#' @details
#' # Getters
#' In addition to what works for `SummarizedExperiment` objects, if `x` is a
#' `TornadoExperiment` object, then:
#' \describe{
#'  \item{`nbin(x)`}{gives the number of bins in `x`. Analogous to `nrow()` and
#'  `ncol()` for the 3rd dimension.}
#'  \item{`binnames(x)`}{gives the names of the bins in `x`. Analogous to
#'    `rownames()` and `colnames()` but for bins}
#'  \item{`binData(x)`}{gets the `binData` slot that contains information about
#'  the 3rd dimension of the `x`.}
#'  \item{`colKey(x); rowKey(x); binKey(x)`}{retrieve the keys associated with
#'  that dimension of `x`.}
#' }
#'
#' # Setters
#' In addition to what works for `SummarizedExperiment` objects, if `x` is a
#' `TornadoExperiment` object, then:
#' \describe{
#'  \item{`binnames(x) <- value`}{sets the names of the bins in `x`.}
#'  \item{`binData(x) <- value`}{sets the `binData` slot if `value` is a
#'  `DataFrame` with `nbin(x)` rows.}
#'  \item{`colKey(x) <- value; rowKey(x) <- value; binKey(x) <- value`}{sets the
#'  keys associated with dimensions of `x` when `value` is a `character(1)`
#'  object and a valid column name in the relevant data slots.}
#' }
#'
#' # Other methods
#' In addition to what works for `SummarizedExperiment` objects, if `x` is a
#' `TornadoExperiment` object, then:
#' \describe{
#'  \item{`x[i, j, k]`}{subsets `x` for the first (`i`), second (`j`) and
#'  third (`k`) dimension.}
#' }
#'
#' @return A `TornadoExperiment` class object.
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' # See ?build_tornado
TornadoExperiment <- function(
  ...,
  binData = DataFrame(),
  colKey = character(),
  rowKey = character(),
  binKey = character()
) {
  sumexp <- SummarizedExperiment(...)
  new_TornadoExperiment(
    sumexp,
    binData = binData,
    colKey = colKey, rowKey = rowKey, binKey = binKey
  )
}

# Assays ------------------------------------------------------------------

# We need the assays method to not be confused with 3 dimnames instead of 2

#' @importFrom SummarizedExperiment assays
setMethod(
  "assays", "TornadoExperiment",
  function(x, withDimnames = TRUE, ...) {
    stopifnot(
      "'withDimnames' must be TRUE or FALSE" =
        isTRUEorFALSE(withDimnames)
    )
    assays <- as(x@assays, "SimpleList")
    if (withDimnames) {
      x_dimnames <- dimnames(x)
      if (is.null(x_dimnames)) {
        x_dimnames <- list(NULL, NULL, NULL)
      }
      assays <- endoapply(
        assays,
        function(a) {
          a_dimnames <- dimnames(a)
          a_dimnames[1:3] <- x_dimnames
          dimnames(a) <- a_dimnames
          a
        }
      )
    }
    assays
  }
)

# binData -----------------------------------------------------------------

## binData ----------------------------------------------------------------

# rowData()/colData() equivalent for bins

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "binData",
  function(x, ...) standardGeneric("binData")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "binData", "TornadoExperiment",
  function(x, ...) {
    x@binData
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "binData<-",
  function(x, ..., value) standardGeneric("binData<-")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setReplaceMethod(
  "binData", "TornadoExperiment",
  function(x, ..., value) {
    if (nrow(value) != nbin(x)) {
      stop("nrow of supplied 'binData' must equal nbin of object.")
    }
    .int$replaceSlots(x, binData = value, check = TRUE)
  }
)

## binnames ---------------------------------------------------------------

# rownames()/colnames() equivalent for bins

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "binnames",
  function(x, ...) standardGeneric("binnames")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "binnames", "TornadoExperiment",
  function(x, ...) {
    rownames(binData(x))
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "binnames<-",
  function(x, ..., value) standardGeneric("binnames<-")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setReplaceMethod(
  "binnames", "TornadoExperiment",
  function(x, value) {
    binData <- binData(x)
    rownames(binData) <- value
    .int$replaceSlots(
      x,
      binData = binData,
      check = FALSE
    )
  }
)

## nbins ------------------------------------------------------------------

# nrow()/ncol() equivalent for bins

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "nbin",
  function(x, ...) standardGeneric("nbin")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "nbin", "TornadoExperiment",
  function(x, ...) {
    dim(x)[3]
  }
)

# Dimensions --------------------------------------------------------------

#' @importFrom SummarizedExperiment colData
#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "dim", "TornadoExperiment",
  function(x) {
    c(length(x), nrow(colData(x)), nrow(binData(x)))
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "dimnames", "TornadoExperiment",
  function(x) {
    ans <- list(names(x), rownames(colData(x)), rownames(binData(x)))

    nulls <- vapply(ans, is.null, logical(1))
    if (all(nulls)) {
      return(NULL)
    }
    ans
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
#' @importFrom SummarizedExperiment colData<-
setReplaceMethod(
  "dimnames", c("TornadoExperiment", "list"),
  function(x, value) {
    dim <- dim(x)
    names <- mapply(
      .int$normarg_names,
      names = value,
      x_class = class(x),
      x_len = dim
    )

    colData <- colData(x)
    rownames(colData(x)) <- names[[2]]

    binData <- binData(x)
    rownames(binData(x)) <- names[[3]]

    .int$replaceSlots(
      x, NAMES = names[[1]], colData = colData, binData = binData,
      check = FALSE
    )
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setReplaceMethod(
  "dimnames", c("TornadoExperiment", "NULL"),
  function(x, value) {
    dimnames(x) <- list(NULL, NULL, NULL)
    x
  }
)

# Subsetting --------------------------------------------------------------

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
#' @importFrom SummarizedExperiment Assays
setMethod(
  "[", "TornadoExperiment",
  function(x, i, j, ..., drop = FALSE) {
    ans <- callNextMethod(x = x, i = i, j = j, drop = drop)

    n_index <- .int$extract_Nindex_from_syscall(sys.call(), parent.frame())

    if (length(n_index) < 3L) {
      return(ans)
    } else {
      k <- n_index[[3]]
    }

    k <- as.integer(NSBS(k, setNames(nbin(x), binnames(x))))

    binData <- binData(x)[k, ]

    assays <- endoapply(assays(ans), function(x) {
      x[, , k, drop = FALSE]
    })
    .int$replaceSlots(
      ans,
      assays  = Assays(assays),
      binData = binData
    )
  }
)

# Keys --------------------------------------------------------------------

## Column keys ------------------------------------------------------------

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "colKey",
  function(x, ...) standardGeneric("colKey")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "colKey", "TornadoExperiment",
  function(x, ...) {
    x@colKey
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "colKey<-",
  function(x, ..., value) standardGeneric("colKey<-")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setReplaceMethod(
  "colKey", c("TornadoExperiment", "character"),
  function(x, value) {
    x@colKey <- value
    validObject(x)
    x
  }
)

## Row keys ---------------------------------------------------------------

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "rowKey",
  function(x, ...) standardGeneric("rowKey")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "rowKey", "TornadoExperiment",
  function(x, ...) {
    x@rowKey
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "rowKey<-",
  function(x, ..., value) standardGeneric("rowKey<-")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setReplaceMethod(
  "rowKey", c("TornadoExperiment", "character"),
  function(x, value) {
    x@rowKey <- value
    validObject(x)
    x
  }
)

## Bin keys ---------------------------------------------------------------

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "binKey",
  function(x, ...) standardGeneric("binKey")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setMethod(
  "binKey", "TornadoExperiment",
  function(x, ...) {
    x@binKey
  }
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setGeneric(
  "binKey<-",
  function(x, ..., value) standardGeneric("binKey<-")
)

#' @export
#' @rdname TornadoExperiment
#' @usage NULL
setReplaceMethod(
  "binKey", c("TornadoExperiment", "character"),
  function(x, value) {
    x@binKey <- value
    validObject(x)
    x
  }
)

# Helpers -----------------------------------------------------------------

get_keys <- function(x, keytype = NULL) {
  if (is.null(keytype)) {
    keytype <- "all"
  }
  switch(
    keytype,
    "row" = rowKey(x),
    "col" = colKey(x),
    "bin" = binKey(x),
    "all" = setNames(c(rowKey(x), colKey(x), binKey(x)),
                     c("row", "col", "bin")),
    stop("Invalid 'keytype'.")
  )
}

get_keydata <- function(x, keytype = "row") {
  switch(
    keytype,
    "row" = rowData(x)[[rowKey(x)]],
    "col" = colData(x)[[colKey(x)]],
    "bin" = binData(x)[[binKey(x)]],
    stop("Invalid 'keytype'.")
  )
}

# Validity ----------------------------------------------------------------

#' @importFrom SummarizedExperiment assay
.valid_TornadoExperiment.assays_bin <- function(object) {
  if (length(object@assays) == 0) {
    return(NULL)
  }
  assays_nbins  <- dim(assay(object, withDimnames = FALSE))[[3]]
  binData_nbins <- nrow(binData(object))
  if (assays_nbins != binData_nbins) {
    txt <- sprintf(
      "\n nb of bins in 'assay' (%d) must equal nb of rows in 'binData' (%d)",
      assays_nbins, binData_nbins
    )
    return(txt)
  }
  NULL
}

#' @importFrom SummarizedExperiment colData
.valid_colkey <- function(object) {
  key <- colKey(object)
  if (length(key) == 0) {
    return(NULL)
  }
  if (!(key %in% colnames(colData(object)))) {
    txt <- paste0(
      "The column key '", key, "' does not exist as a column name in ",
      "the 'colData'."
    )
    return(txt)
  }
  NULL
}

#' @importFrom SummarizedExperiment rowData
.valid_rowkey <- function(object) {
  key <- rowKey(object)
  if (length(key) == 0) {
    return(NULL)
  }
  if (!(key %in% colnames(rowData(object)))) {
    txt <- paste0(
      "The row key '", key, "' does not exist as a column name in ",
      "the 'rowData'."
    )
    return(txt)
  }
  NULL
}

.valid_binkey <- function(object) {
  key <- binKey(object)
  if (length(key) == 0) {
    return(NULL)
  }
  if (!(key %in% colnames(binData(object)))) {
    txt <- paste0(
      "The bin key '", key, "' does not exist as a column name in ",
      "the 'binData'."
    )
    return(txt)
  }
  NULL
}

setValidity2(
  "TornadoExperiment",
  function(object) {
    msg <- NULL
    msg <- c(msg, .valid_TornadoExperiment.assays_bin(object))
    msg <- c(msg, .valid_colkey(object))
    msg <- c(msg, .valid_rowkey(object))
    msg <- c(msg, .valid_binkey(object))
    if (length(msg)) {
      return(msg)
    }
    TRUE
  }
)

# Show --------------------------------------------------------------------

setMethod(
  "show", "TornadoExperiment",
  function(object) {
    callNextMethod()

    binnames <- binnames(object)
    if (!is.null(binnames)) {
      coolcat("binnames(%d): %s\n", binnames)
    } else {
      cat("binnames: NULL\n")
    }
    coolcat("binData names(%d): %s\n", names(binData(object)))
  }
)
