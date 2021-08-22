# Formatting --------------------------------------------------------------

# These are a bunch of functions that help format data a TornadoExperiment in
# the `build_tornado` function. These are all intended for internal use and
# aren't exported.

# Sample data -------------------------------------------------------------

format_sample_data <- function(data, data_arg, n, ..., barcode_groups = NULL) {
  file_classes <- c("RTLFile", "RTLFileList",
                    "RsamtoolsFile", "RsamtoolsFileList",
                    "BiocFile", "BiocFileList")
  # Set names
  nms <- format_sample_name(data, data_arg, n, names = NULL,
                            barcode_groups = barcode_groups, ...)

  if (hasMethod("path", class(data))) {
    files <- path(data)
  } else if (is.character(data) && all(file.exists(data))) {
    files <- data
  } else {
    files <- NULL
  }

  if (!is.null(barcode_groups)) {
    barcode_groups <- as(barcode_groups, "CharacterList")
    if (all(names(barcode_groups) %in% nms)) {
      barcode_groups <- barcode_groups[nms]
    }
  }

  df_args <- list(
    sample_name = nms,
    argument = data_arg,
    file = files,
    barcodes = barcode_groups,
    row.names = nms
  )
  df_args <- df_args[lengths(df_args) > 0]

  do.call(DataFrame, df_args)
}

# Feature data ------------------------------------------------------------

format_feature_data <- function(features, feat_arg, n) {
  lens <- metadata(features)$lengths
  metadata(features)$lengths <- NULL

  set <- lens %||% Rle(1, length(features))

  df_args <- list(
    set = set,
    argument = Rle(feat_arg, length(set))
  )
  df_args <- df_args[lengths(df_args) > 0]

  mcols(features) <- do.call(DataFrame, df_args)
  features
}

# Bin data ----------------------------------------------------------------

#' @importFrom IRanges shift IRanges
#' @importFrom BiocGenerics start end
format_bin_data <- function(bins) {
  len <- length(bins)
  DataFrame(
    bin_id = runValue(bins),
    range  = shift(IRanges(start(bins), end(bins)), shift = - 0.5 * len)
  )
}

# Helpers -----------------------------------------------------------------

# TODO: Include BiocFile when switching to BioC 1.13
#' @importClassesFrom rtracklayer RTLFileList
#' @importClassesFrom Rsamtools RsamtoolsFileList
setClassUnion(
  "FileListUnion",
  c("RTLFileList", "RsamtoolsFileList")
)

#' @importClassesFrom rtracklayer RTLFile
#' @importClassesFrom Rsamtools RsamtoolsFile
setClassUnion(
  "SingleFileUnion",
  c("RTLFile", "RsamtoolsFile")
)

# Sample names ------------------------------------------------------------

setGeneric(
  "format_sample_name",
  function(data, data_arg, n, names = NULL, ...)
    standardGeneric("format_sample_name"),
  signature = "data"
)

setMethod(
  "format_sample_name", "TabixFile",
  function(data, data_arg, n, names = NULL, barcode_groups = NULL, ...) {
    names <- names %||% names(barcode_groups)
    callNextMethod(data = data, data_arg = data_arg, n = n, names = names, ...)
  }
)

#' @importFrom BiocGenerics path
setMethod(
  "format_sample_name", "FileListUnion",
  function(data, data_arg, n, names = NULL, ...) {
    names <- names %||% names(data) %||% format_filename(path(data))
    callNextMethod(data = data, data_arg = data_arg, n = n, names = names, ...)
  }
)

setMethod(
  "format_sample_name", "SingleFileUnion",
  function(data, data_arg, n, names = NULL, ...) {
    names <- names %||% format_filename(path(data))
    callNextMethod(data = data, data_arg = data_arg, n = n, names = names, ...)
  }
)

setMethod(
  "format_sample_name", "character",
  function(data, data_arg, n, names = NULL, ...) {
    names <- names %||% names(data) %||% format_filename(data)
    callNextMethod(data = data, data_arg = data_arg, n = n, names = names, ...)
  }
)

setMethod(
  "format_sample_name", "ANY",
  function(data, data_arg, n, names = NULL, ...) {
    if (hasMethod("path", class(data))) {
      if (length(data) == 1) {
        names <- names %||% format_filename(path(data))
      } else {
        names <- names %||% names(data) %||% format_filename(path(data))
      }
    }
    names <- names %||% data_arg
    if (length(names) != n) {
      names <- paste(names, seq_len(n), sep = "_")
    }
    names
  }
)

format_filename <- function(path) {
  # Regex from: https://stackoverflow.com/a/63939003/11374827
  gsub("\\..[^\\.]*$", "", basename(path))
}

