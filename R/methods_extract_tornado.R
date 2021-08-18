# Generic -----------------------------------------------------------------

setGeneric(
  "extract_tornado",
  signature = "data",
  function(
    data, features, bins,
    pad_value = 0, ...
  ) standardGeneric("extract_tornado")
)

# Bigwig Files ------------------------------------------------------------

# Abuse magic summary method for bigwig files
#' @importClassesFrom rtracklayer BigWigFileList
setMethod(
  "extract_tornado",
  signature = c(data = "BigWigFileList"),
  function(
    data, features, bins,
    pad_value = 0, ...
  ) {
    ans <- matrix(0, nrow = length(features), ncol = nrun(bins))
    ans <- vapply(
      data, summary, ans,
      which = features, size = nrun(bins),
      as = "matrix", type = "mean",
      defaultValue = pad_value
    )
    aperm(ans, c(2, 1, 3))
  }
)

# Convert BigWigFile to BigWigFileList
#' @importClassesFrom rtracklayer BigWigFile
setMethod(
  "extract_tornado",
  signature = c(data = "BigWigFile"),
  function(
    data, features, bins,
    pad_value = 0, ...
  ) {
    nms <- names(data)
    data <- BigWigFileList(setNames(path(data), nms))
    callGeneric(data, features, bins, pad_value, ...)
  }
)

# Try to coerce characters to BigWigFileLists
setMethod(
  "extract_tornado",
  signature = c(data = "character"),
  function(
    data, features, bins,
    pad_value = 0, ...
  ) {
    if (!all(file.exists(data))) {
      stop("The `data` argument doesn't appear to be a valid path to a file.",
           call. = FALSE)
    }
    if (!all(endsWith(data, ".bw"))) {
      stop("The file type is unsupported.",
           call. = FALSE)
    }
    data <- BigWigFileList(data)
    callGeneric(data, features, bins, pad_value, ...)
  }
)

# Genomic Ranges ----------------------------------------------------------

setOldClass("ListOfRleList")

#' @importClassesFrom GenomicRanges GRanges GNCList
#' @importFrom GenomicRanges coverage GNCList
#' @importFrom IRanges overlapsAny
setMethod(
  "extract_tornado",
  signature = c(data = "GRanges"),
  function(
    data, features, bins,
    pad_value = 0, ...
  ) {
    data <- data[overlapsAny(data, GNCList(reduce(features)))]
    data <- coverage(data)
    data <- structure(list(data), class = "ListOfRleList")
    callGeneric(data, features, bins, pad_value, ...)
  }
)

#' @importClassesFrom GenomicRanges GRangesList
setMethod(
  "extract_tornado",
  signature = c(data = "GRangesList"),
  function(
    data, features, bins,
    pad_value = 0, ...
  ) {
    data <- stack(data)
    data <- data[overlapsAny(data, GNCList(reduce(features)))]
    data <- lapply(split(data, data$name), coverage)
    data <- structure(data, class = "ListOfRleList")
    callGeneric(data, features, bins, pad_value, ...)
  }
)

# TabixFile ---------------------------------------------------------------

#' @importClassesFrom Rsamtools TabixFile
#' @importFrom Rsamtools headerTabix
setMethod(
  "extract_tornado",
  signature = c(data = "TabixFile"),
  function(
    data, features, bins,
    pad_value = 0,
    barcode_column = NULL,
    barcode_groups = NULL,
    ...
  ) {
    if (is.null(barcode_column) || is.null(barcode_groups)) {
      barcode_column <- barcode_groups <- NULL
    }
    if (!is.null(barcode_groups)) {
      groups <- names(barcode_groups) %||% seq_along(barcode_groups)
      groups <- rep.int(groups, lengths(barcode_groups))
      barcodes <- unlist(barcode_groups, use.names = FALSE)
    } else {
      barcodes <- NULL
    }

    # Setup column classes from tabix header
    info <- headerTabix(data)$indexColumns
    colclasses <- rep(NA, max(info, barcode_column))
    colclasses[info] <- c("character", "integer", "integer")
    colclasses[barcode_column] <- "character"

    # Read in Tabix data
    data <- unlist(Rsamtools::scanTabix(data, param = features),
                   use.names = FALSE)
    data <- read.table(
      text = data, sep = "\t", quote = "", comment.char = "",
      colClasses = colclasses, check.names = FALSE, encoding = "UTF8",
      as.is = TRUE
    )

    # Generate index for splitting
    if (!is.null(barcodes)) {
      idx <- groups[match(data[[barcode_column]], barcodes)]
    } else {
      idx <- rep(1L, nrow(data))
    }

    data <- GRanges(
      data[[info[["seq"]]]],
      IRanges(
        data[[info[["start"]]]],
        data[[info[["end"]]]]
      )
    )

    data <- structure(lapply(split(data, idx), coverage),
                      class = "ListOfRleList")
    callGeneric(data, features, bins, pad_value, ...)
  }
)

# Rle lists ---------------------------------------------------------------

# Essentially does all the non-bigwig file operations

setMethod(
  "extract_tornado",
  signature = c(data = "ListOfRleList"),
  function(
    data, features, bins,
    pad_value = 0, ...
  ) {
    data <- pad_coverage(data, features, pad_value)

    # Calculate indices
    feat_len <- length(features)
    index <- rep(end(bins), feat_len)
    index <- index + rep(0:(feat_len - 1) * length(bins), each = nrun(bins))

    # Summarise coverage
    template <- matrix(0, nrow = nrun(bins), ncol = length(features))
    ans <- vapply(data, summarise_core, template,
                  feats = features, index = index, dim = dim(template))
    ans / runLength(bins)[1]
  }
)



# Helpers -----------------------------------------------------------------

# Old wrapper, replaced by ListOfRleList objects
# summarise_coverage <- function(coverage, feats, bins) {
#   template <- matrix(0, nrow = nrun(bins), ncol = length(feats))
#
#   # Calculate indices
#   index <- rep(end(bins), length(feats))
#   index <- index + rep(0:(length(feats) - 1) * length(bins), each = nrun(bins))
#
#   ans <- vapply(coverage, summarise_workhorse, template,
#                 feats = feats, index = index, dim = dim(template))
#
#   # Flip minus strand orientation
#   swap <- which(decode(strand(feats) == "-"))
#   ans[ , swap, ] <- ans[rev(seq_len(nrun(bins))), swap, ]
#
#   ans / runLength(bins)[1]
# }

#' @importFrom GenomeInfoDb seqnames seqlevelsStyle
#' @importMethodsFrom GenomicRanges range
pad_coverage <- function(coverage, features, pad_value = 0) {
  # Get 'empirical' seqlengths from features
  features <- range(features, ignore.strand = TRUE)
  features <- setNames(end(features), decode(seqnames(features)))

  covnames <- Reduce(intersect, lapply(coverage, names))
  padme    <- check_seqlevels(covnames, unique(names(features)))

  lapply(coverage, function(covr) {
    seqs <- intersect(names(covr), names(features))
    pad  <- features[seqs] - lengths(covr)[seqs]
    pad  <- pad[seqs]
    for (i in names(pad)[pad > 0]) {
      covr[[i]] <- c(covr[[i]], Rle(pad_value, pad[i]))
    }
    for (i in names(padme)) {
      covr[[i]] <- Rle(pad_value, features[i])
    }
    return(covr)
  })
}

check_seqlevels <- function(seqlevels_x, seqlevels_y) {
  if (length(intersect(seqlevels_x, seqlevels_y)) == 0) {
    style_cover <- seqlevelsStyle(seqlevels_x)
    style_feats <- seqlevelsStyle(seqlevels_y)
    if (style_feats != style_cover) {
      msg <- paste0(
        "The features and data have no common sequences. This might be due ",
        "to a difference in style: the features follow the '", style_feats,
        "' convention, whereas the data follow the '", style_cover,
        "' convention. Try adjusting the features with ",
        "`GenomeInfoDb::seqlevelsStyle()`."
      )
    } else {
      msg <- paste0(
        "The features and data have no common sequences, even though they ",
        "share style conventions of naming sequences."
      )
    }
    stop(msg, call. = FALSE)
  }
  dif <- setdiff(seqlevels_y, seqlevels_x)
  if (length(dif)) {
    if (length(dif) == 1) {
      difnames <- dif
    } else {
      difnames <- paste0(
        paste0(dif[-length(dif)], collapse = ", "),
        " and ", dif[length(dif)]
      )
    }
    msg <- paste0(
      "The following sequence names were missing from the data but present in ",
      "the features: ", difnames, ". The data will be padded for these ",
      "sequences."
    )
  }
  dif
}

# Core function -----------------------------------------------------------

# This function exists for cases when the input is an RleList object
#
# There is a trick in `summarise_coverage` inspired by the 'summed area table'
# but in just 1 dimension (https://en.wikipedia.org/wiki/Summed-area_table).
# Instead of calculating a bunch of means over bins, the cumulative sum is
# calculated, which is later split into means.
#
# Example:
# # Classic split along bins, calculate sums, divide by lengths
# classic <- function(x, bins) {
#   unname(vapply(split(x, bins), sum, numeric(1))) / runLength(bins)[1]
# }
#
# # Only consider ends of cumulative sum
# trick <- function(x, bins) {
#   x <- cumsum(x)[end(bins)]
#   x[-1] <- x[-1] - x[-length(x)]
#   x / runLength(bins)[1]
# }
#
# # Test data
# x <- rnorm(1e6)
# bins <- Rle(1:1e4, 100)
#
# # Compare
# bench::mark(
#   {classic(x, bins)},
#   {trick(x, bins)},
# )
# summarise_core <- function(coverage, feats, index, dim) {
#   # Select coverage at features
#   coverage <- coverage[feats]
#
#   # Calculate cumulative sums
#   coverage <- unlist(coverage, use.names = FALSE, recursive = FALSE)
#   coverage <- cumsum(rep.int(runValue(coverage), runLength(coverage)))
#
#   # Take bin ends and subtract previous bin end
#   coverage <- coverage[index]
#   coverage <- coverage - c(0, coverage[-length(coverage)])
#
#   # Reshape as matrix
#   dim(coverage) <- dim
#   coverage
# }


#' Title
#'
#' @param data A named `RleList` object, wherein names correspond to
#'   `seqlevels(feats)`.
#' @param feats A `GRanges` object of equal width wherein the `seqlevels()`
#'   correspond to the names in `data`.
#' @param index An `integer` vector with bin ends.
#' @param dim An `integer(2)` with the output dimensions.
#'
#' @return A `matrix` with dimensions
#' @noRd
summarise_core <- function(data, feats, index, dim) {
  # Select coverage at features
  data <- data[feats]
  data <- unlist(data, FALSE, FALSE)

  # Find which run index falls into and how many elements into that run
  coverend <- end(data)
  interval <- findInterval(index, coverend, left.open = TRUE) + 1
  overflow <- index - coverend[interval]

  # Calculate cumsums
  ans <- cumsum(runValue(data) * runLength(data))[interval]
  ans <- ans + overflow * runValue(data)[interval]

  # Convert cumsums to binned averages
  ans <- ans - c(0, ans[-length(ans)])
  dim(ans) <- dim
  ans
}

