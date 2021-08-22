#' Generation of test data
#'
#' For the purposes of examples and testing the package, tornadoplot includes
#' some functions that make up some dummy data that isn't too large.
#'
#' @param pattern A `character(1)` giving an initial part of a file name.
#' @param dir A `character(1)` giving a directory name where to create the test
#'   file.
#'
#' @details The `dummy_granges_data()` and `dummy_bigwig()` generators make
#'   their data around the genomic ranges returned by `dummy_features()`.
#'
#'   The `dummy_bigwig()` and `dummy_tabix()` files create files on your system.
#'   As their usage is intended as temporary files, don't forget to `unlink()`
#'   their paths and indices after they are no longer needed.
#'
#'   The `dummy_tornado()` makes a relatively small
#'   \linkS4class{TornadoExperiment}.
#'
#' @return The return values differ per function. \describe{
#'  \item{`dummy_features`}{A `GRanges` object of length 2.}
#'  \item{`dummy_granges_data`}{A `GRanges` object of length 20.}
#'  \item{`dummy_bigwig`}{A `BigWigFile` object.}
#'  \item{`dummy_tabix`}{A `TabixFile` object with `"A"` and `"B"` in the 4th
#'  column.}
#'  \item{`dummy_tornado`}{A \linkS4class{TornadoExperiment} object with 2
#'  samples, 4 features and 50 bins.}
#' }
#'
#' @name testdata
#'
#' @examples
#' # All the functions run without arguments
#' feats <- dummy_features()
#' gr    <- dummy_granges_data()
#' bw    <- dummy_bigwig()
#' tbx   <- dummy_tabix()
#'
#' # Don't forget to clean up the temporary files when you're done
#' unlink(c(path(bw), path(tbx), Rsamtools::index(tbx)))
NULL

#' @rdname testdata
#' @export
#' @importFrom GenomicRanges GRanges
dummy_features <- function() {
  GRanges(c("chr2:501-2501", "chr9:10501-12501"))
}

#' @rdname testdata
#' @importFrom GenomeInfoDb Seqinfo
#' @export
dummy_granges_data <- function() {
  GRanges(
    rep(c("chr2", "chr9"), each = 10),
    IRanges(
      c(seq(1, 1450, by = 161),
        seq(10001, 11450, by = 161)),
      width = rep(seq(2998, 100, by = -322), 2)
    ),
    seqinfo = Seqinfo(c("chr2", "chr9"), c(182113224, 124595110))
  )
}

#' @rdname testdata
#' @export
#' @importFrom Rsamtools bgzip indexTabix TabixFileList
dummy_tabix <- function(pattern = "file", dir = tempdir()) {
  file <- tempfile(pattern, tmpdir = dir, fileext = ".bed")
  on.exit(unlink(file))
  gr <- dummy_granges_data()
  cells <- rep_len(c("A", "B"), length(gr))
  txt <- paste(
    as.character(seqnames(gr)),
    start(gr), end(gr), cells,
    sep = "\t"
  )
  cat(txt, file = file, sep = "\n")
  newfile <- Rsamtools::bgzip(file, paste0(file, ".gz"))
  idx <- Rsamtools::indexTabix(newfile, seq = 1, start = 2, end = 3)
  Rsamtools::TabixFile(setNames(newfile, "test tabix"), idx)
}

#' @rdname testdata
#' @export
#' @importFrom rtracklayer export.bw BigWigFile
dummy_bigwig <- function(pattern = "file", dir = tempdir()) {
  file <- tempfile(pattern, tmpdir = dir, fileext = ".bw")
  gr <- GRanges(
    rep(c("chr2", "chr9"), each = 10),
    IRanges(
      c(seq(1, 2701, by = 300),
        seq(10001, 12701, by = 300)),
      width = 300
    ),
    score = c(dlaplace(1:10, 4, 1),
              dlaplace(1:10, 6, 2)),
    seqinfo = Seqinfo(c("chr2", "chr9"), c(182113224, 124595110))
  )
  rtracklayer::export.bw(gr, file)
  BigWigFile(setNames(file, "test bigwig"))
}

#' @rdname testdata
#' @export
dummy_tornado <- function() {
  b <- resolve_bins(2000, nbin = 50)
  f <- GRanges(c("chr1:1001-3000", "chr2:1001-3000",
                 "chr3:1001-3000", "chr4:1001-3000"))
  f <- split(f, c("A", "A", "B", "B"))
  f <- resolve_features(f, b)

  d <- c("dummy_ctrl.bw", "dummy_treat")

  ans <- array(NA_real_, dim = c(length(f), length(d), nrun(b)))
  ans[] <- dlaplace(
    as.vector(slice.index(ans, 3)),
    mu = (nrun(b) + 1)/2, b = as.vector(slice.index(ans, 1)) * 10
  ) * as.vector(slice.index(ans, 2))

  TornadoExperiment(
    assays = SimpleList(tornado = ans),
    rowRanges = format_feature_data(f, "feats", n = dim(ans)[1]),
    colData   = format_sample_data(d, "files", n = dim(ans)[2]),
    binData   = format_bin_data(b),
    rowKey = "set", colKey = "sample_name", binKey = "bin_id"
  )
}

#' Sox2 and Oct4 binding sites
#'
#' This dataset contains a sample of genomic binding sites of the Sox2 and Oct4
#' proteins in mouse embryonic stem cells.
#'
#' @usage NULL
#' @format A \linkS4class{GRangesList} with 3 elements, each having 500 genomic
#' locations of binging sites.
#' \describe{
#'   \item{Sox2}{500 genomic loci where Sox2 binds}
#'   \item{Oct4}{500 genomic loci where Oct4 binds}
#'   \item{Both}{500 genomic loci where both Sox2 and Oct4 bind}
#' }
#'
#' @details The genomic loci were obtained by chromatin immunoprecipitation
#' followed by sequencing (ChIP-seq) with antibodies against the Sox2 and Oct4
#' proteins. The model system used were wildtype mouse embryonic stem cells.
#'
#' The original data was deposited in the gene expression omnibus (GEO) under
#' accession
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74112}{GSE74112}.
#'
#' The data were aligned against the "mm10" reference genome by Cistrome. Other
#' preprocessing steps are described at the
#' \href{http://cistrome.org/db/#/about}{Cistrome website}.
#'
#' Finally, locations in standard chromosomes were reduced to unique sites and
#' categorised based on overlap with Sox2 binding sites and Oct4 binding sites
#' or both. From each of the three resultant categories, 500 locations were
#' randomly sampled without replacement.
#'
#' @references
#' Liu, Ziying *et al*. Catalytic-Independent Functions of PARP-1 Determine Sox2
#' Pioneer Activity at Intractable Genomic Loci (2017) Molecular Cell *65* **4**
#' p589-603 https://doi.org/10.1016/j.molcel.2017.01.017
#'
#' @source
#' The files used to construct this data were downloaded as bed files from
#' \url{http://cistrome.org/db/#/} (accessed at 2021-08-18). The files can be
#' unambiguously found by entering the following numbers in the search bar:
#' "73468" for the Sox2 sample and "73466" for the Oct4 sample.
"octsox_peaks"
