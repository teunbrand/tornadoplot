# Libraries ---------------------------------------------------------------

library(GenomicRanges)
library(rtracklayer)

# Files -------------------------------------------------------------------

# These files were downloaded from cistrome: http://cistrome.org/db/#/
# They are bed files from data deposited in GSE74112 from PMID: 28212747
# The data is chromatin immunoprecipitation followed by sequencing (ChIP-seq),
# where the proteins Sox2 and Oct4/Pou5f1 have been precipitated from
# mouse embryonic stem cells (mESC).

dir <- file.path("/DATA", "users", "t.vd.brand", "test_data")
files <- file.path(
  dir,
  c("73466_peaks_oct4.bed",
    "73468_peaks_sox2.bed")
)

# Data import -------------------------------------------------------------

si <- SeqinfoForUCSCGenome("mm10")
si <- keepStandardChromosomes(si, "Mus_musculus")

# We don't want any of the scaffold sequences

bed <- lapply(files, import)
bed <- lapply(bed, keepStandardChromosomes, "Mus_musculus", "coarse")
bed <- as(bed, "GRangesList")
bed <- setNames(bed, c("oct4", "sox2"))

# We want to find sites where Sox2 and Oct4 co-bind and where they bind
# uniquely.

red <- reduce(stack(bed))
red$sox2 <- overlapsAny(red, bed$sox2)
red$oct4 <- overlapsAny(red, bed$oct4)
red$cat <- ifelse(red$sox2 & red$oct4, "both",
                  ifelse(red$sox2, "sox2", "oct4"))
red <- GRanges(seqnames(red), IRanges(start(red), end(red)),
               seqinfo = si, cat = red$cat)

# Choosing sites ----------------------------------------------------------

set.seed(0)
octsox_peaks <- GRangesList(
  Sox2 = sample(red[red$cat == "sox2"], 500),
  Oct4 = sample(red[red$cat == "oct4"], 500),
  Both = sample(red[red$cat == "both"], 500)
)
octsox_peaks <- sort(octsox_peaks)

usethis::use_data(octsox_peaks, overwrite = TRUE)
