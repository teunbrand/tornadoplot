setGeneric(
  "extract_coverage", signature = "data",
  function(data, selection) standardGeneric("extract_coverage")
)

#' @importClassesFrom GenomicRanges GRanges
setMethod(
  "extract_coverage", signature = c(data = "GRanges"),
  function(data, selection) {
    data <- data[overlapsAny(data, GNCList(selection))]
    list(coverage(data))
  }
)

#' @importClassesFrom GenomicRanges GRangesList
setMethod(
  "extract_coverage", signature = c(data = "GRangesList"),
  function(data, selection) {
    data <- stack(data)
    data <- data[overlapsAny(data, GNCList(selection))]
    lapply(split(data, data$name), coverage)
  }
)

#' @importClassesFrom rtracklayer BigWigFile
#' @importFrom rtracklayer import.bw BigWigSelection
setMethod(
  "extract_coverage", signature = c(data = "BigWigFile"),
  function(data, selection) {
    selection <- BigWigSelection(selection)
    list(import.bw(data, selection = selection, as = "RleList"))
  }
)

#' @importClassesFrom rtracklayer BigWigFileList
#' @importFrom rtracklayer import.bw
setMethod(
  "extract_coverage", signature = c(data = "BigWigFileList"),
  function(data, selection) {
    lapply(data, import.bw, selection = selection, which = selection, as = "RleList")
  }
)

#' @importFrom rtracklayer BigWigFileList
setMethod(
  "extract_coverage", signature = c(data = "character"),
  function(data, selection) {
    if (!all(file.exists(data))) {
      stop("The `data` argument doesn't appear to be a valid path to a file.",
           call. = FALSE)
    }
    if (!all(endsWith(data, ".bw"))) {
      stop("The file type is unsupported.",
           call. = FALSE)
    }
    data <- BigWigFileList(data)
    extract_coverage(data, selection)
  }
)
