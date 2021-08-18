#' Acquire tornado plot data
#'
#' Aggregates coverage data in an array. Coverage can be read in from
#' \linkS4class{BigWigFile}s or computed from a \linkS4class{GRanges} object or
#' a \linkS4class{TabixFile}.
#'
#' @param features A \linkS4class{GRanges} or \linkS4class{GRangesList} object
#'   containing genomic loci of interest.
#' @param data The data source to get coverages. Either a \linkS4class{GRanges},
#'   \linkS4class{GRangesList}, \linkS4class{TabixFile},
#'   \linkS4class{BigWigFile}, \linkS4class{BigWigFileList} object or
#'   `character` path to `.bw` files.
#' @param width,binwidth,nbin An `integer(1)` in basepairs with respectively a
#'   common feature width to centre features in, a size of bins to summarise
#'   coverage in or the number of bins to summarise coverage in. Only two of
#'   these need to be defined.
#' @param ... Arguments used for \linkS4class{TabixFile} input (see single cell
#' data section in details):
#' \describe{
#'  \item{`barcode_column`}{An `integer(1)` for which column in the
#'  \linkS4class{TabixFile} are barcodes, typically describing cells in single
#'  cell assays.}
#'  \item{`barcode_groups`}{A (named) `list` of barcodes wherein every element
#'  is a `character` vector with barcodes that belong to the same groups,
#'  typically clusters of cells.}
#' }
#' @param pad_value A `numeric(1)` to use for padding when
#'   [`seqlenghts(features)`][GenomeInfoDb::seqlengths()] are unknown or
#'   features exceed known sequence lengths.

#' @details
#' ## Features
#' The genomic ranges provided as the `features` argument get
#' [`resize`][GenomicRanges::intra-range-methods]d to all have widths equal to
#' the (computed) `width` argument, with fixed centres. Features with a negative
#' strand will have their output flipped in the output, such that stranded
#' features yield 5' -> 3' coverage data. If this is undesired, the features can
#' be [`unstrand`][BiocGenerics::unstrand()]ed to prevent this flipping. When
#' the `features` are provided as a \linkS4class{GRangesList}, every list
#' element is interpreted as part of a different feature set. Internally these
#' get unlisted an their set membership is tracked as the `feature_set` column
#' in the `rowData` slot.
#'
#' ## Data
#'
#' ### Bulk data
#' If the data is a (set of) bigwig files, the output is constructed through the
#' [`summary()`][rtracklayer::BigWigFile-class] method for bigwig files. For
#' other types of data, the [`coverage()`][GenomicRanges::coverage-methods] is
#' computed at the locations of the features. This basepair resolution coverage
#' is subsequently binned and averaged.
#'
#' ### Single cell data
#' A popular format for storing single cell chromatin data is as a
#' \linkS4class{TabixFile}. For example, the 10x Genomics 'cellranger' pipeline
#' for single cell ATAC-seq produces a
#' [fragments](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments)
#' tabix file, wherein the 4th column indicates the barcode of a
#' cell. When the `data` argument is a tabix file, the `barcode_column` argument
#' instructs this function where to look for barcodes. The list elements of the
#' `barcode_groups` argument can then be matched against that column in the
#' tabix file to extract the data for every group of cells. Next, the coverage
#' for every group of cells is calculated, binned and averaged in the same way
#' bulk data is.

#' @return A \linkS4class{TornadoExperiment} object with the following populated
#' slots:
#' \describe{
#'   \item{assays}{Has an `n features` &times; `m samples` &times; `o bins` 3D
#'   array with coverage data.}
#'   \item{colData}{Information that could be derived from the `data` argument.}
#'   \item{rowRanges}{Flattened and resized `GRanges` derived from the
#'   `features` argument}
#'   \item{rowData}{Has a `feature_set` column. See the 'Features' details
#'   subsection.}
#'   \item{binData}{A \linkS4class{DataFrame} with a `bin_id` and `range`
#'   column.}
#' }
#' @export
#' @importFrom BiocGenerics strand
#'
#' @examples
#' # Some very small features and data that works
#' feats <- dummy_features()
#' dat   <- dummy_granges_data()
#'
#' # Make a tornado
#' tor   <- build_tornado(feats, dat, width = 2000)
#'
#' # Plotting the tornado
#' autoplot(tor)
build_tornado <- function(
  features, data,
  width = 4000, binwidth = 25, nbin = NULL,
  ...,
  pad_value = 0
) {
  # Capture argument names
  data_arg <- deparse(substitute(data))
  feat_arg <- deparse(substitute(features))

  # Resolve arguments
  bins     <- resolve_bins(width, binwidth, nbin)
  features <- resolve_features(features, bins)

  # Extract coverage from data
  ans <- extract_tornado(data, features, bins, pad_value = 0, ...)

  # Flip minus strand orientation
  swap <- which(decode(strand(features) == "-"))
  if (length(swap)) {
    ans[ , swap, ] <- ans[rev(seq_len(nrun(bins))), swap, ]
  }
  ans <- aperm(ans, c(2, 3, 1))

  # Format answer
  TornadoExperiment(
    assays    = SimpleList(tornado = ans),
    rowRanges = format_feature_data(features, feat_arg, n = dim(ans)[1]),
    colData   = format_sample_data(data, data_arg, n = dim(ans)[2], ...),
    binData   = format_bin_data(bins),
    rowKey = "set", colKey = "sample_name", binKey = "bin_id"
  )
}

# Resolve arguments -------------------------------------------------------

#' @importFrom S4Vectors Rle
resolve_bins <- function(width = NULL, binwidth = NULL, nbin = NULL) {
  # Check for null and coerce to integer
  width    <- norm_type_or_null(width, integer())
  binwidth <- norm_type_or_null(binwidth, integer())
  width    <- norm_type_or_null(width, integer())

  check_null <- c(is.null(width), is.null(binwidth), is.null(nbin))
  if (sum(check_null) > 1) {
    stop("Define 2 out of 3: `width`, `binwidth` and `nbin`.",
         call. = FALSE)
  }
  if (sum(check_null) == 0) {
    # A defined nbin overrules a defined binwidth
    binwidth <- width / nbin
  }

  if (check_null[3]) {
    nbin <- width / binwidth
    if (nbin %% 1 != 0) {
      stop("`width` is not a multiple of binwidth.")
    }
  }
  if (check_null[1]) {
    width <- nbin * binwidth
  }

  if (check_null[2]) {
    binwidth <- width / nbin
    if (binwidth %% 1 != 0) {
      stop("`width` does not fit an integer of `nbin`.",
           call. = FALSE)
    }
  }

  Rle(seq_len(nbin), binwidth)
}

#' @importFrom GenomicRanges resize granges reduce
#' @importFrom GenomeInfoDb keepSeqlevels seqlevelsInUse
resolve_features <- function(features, bins,
                             arg = deparse(substitute(features))) {
  # Flatten GRangesList and note lengths
  if (inherits(features, "GRangesList")) {
    feat_len <- names(features) %||% seq_along(features)
    feat_len <- Rle(feat_len, lengths(features))
    features <- unlist(features)
  } else {
    feat_len <- NULL
  }
  # At this point we should have a GRanges object
  if (!inherits(features, "GRanges")) {
    msg <- paste0(
      "Expected the `", arg,
      "` argument be a `GRanges` or `GRangesList` object."
    )
    stop(msg, .call = FALSE)
  }

  # Set common width for all features
  features <- resize(granges(features), length(bins), fix = "center")
  # Drop unused seqlevels
  features <- keepSeqlevels(features, seqlevelsInUse(features))

  if (length(features) > 50000) {
    # Might be a ridiculous request
    flen <- length(features)
    pages <- round(((flen / 300) * 2.54) / 29.7)
    warning(paste0(
      "Printing ", flen, " features at 1 pixel-width each, would require about ",
      pages, " full A4 sized papers, at 300 ppi resolution.\nIf building the ",
      "tornado takes a long time, consider downsampling the features."
    ))
  }
  # Attach feature length metadata
  metadata(features) <- list(lengths = feat_len)
  return(features)
}

resolve_selection <- function(features, arg = deparse(substitute(features))) {
  if (inherits(features, "GRangesList")) {
    features <- unlist(features)
  }
  if (!inherits(features, "GRanges")) {
    msg <- paste0("Expected `", arg,
                  "` to be a `GRanges` or `GRangesList` object.")
    stop(msg, .call = FALSE)
  }
  reduce(features)
  features <- keepSeqlevels(features, seqlevelsInUse(features))
}






