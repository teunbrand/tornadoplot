% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shape_tornado.R
\name{sort_tornado}
\alias{sort_tornado}
\title{Reorder features}
\usage{
sort_tornado(
  tornado,
  assay_name = assayNames(tornado),
  bin_subset = NULL,
  bin_weights = NULL,
  sample_subset = NULL,
  decreasing = TRUE,
  ...
)
}
\arguments{
\item{tornado}{A \linkS4class{TornadoExperiment} object.}

\item{assay_name}{A \code{character(1)}: one of the assay names.}

\item{bin_subset, bin_weights}{Define either \code{bin_subset} or \code{bin_weights} (or
neither), but not both. \code{bin_subset} should be a \code{logical} or \code{numeric}
vector with subsetting indices for the number of bins. \code{bin_weights} should
be a \code{numeric} vector of weights, the same length as there are bins. If
both \code{bin_subset} and \code{bin_weights} are \code{NULL} (default), the weights are
constructed with the density function of the Laplace distribution wherein
the location is centred at the middle of the bins and the scale parameter
is 10\% of the number of bins.}

\item{sample_subset}{A \code{character}, \code{logical} or \code{integer} vector for
subsetting the samples.}

\item{decreasing}{A \code{logical(1)}: sort from high to low (\code{TRUE}) or from low
to high (\code{FALSE})?}

\item{...}{Not currently in use.}
}
\value{
The \code{tornado} dataset with the features reordered.
}
\description{
The \code{sort_tornado()} function reorders the features of a
\linkS4class{TornadoExperiment} object based on the weighted sum across all
bins.
}
\examples{
# A tornado that isn't sorted
tor <- dummy_tornado()[c(2,1,4,3),]

# Sorting based on a particular sample
x <- sort_tornado(tor, sample_subset = "dummy_ctrl")

# Sorting based on a subset of bins
x <- sort_tornado(tor, bin_subset = 21:30)

# Sorting with custom weights
nbin <- nbin(tor)
w <- dnorm(seq_len(nbin), mean = (nbin + 1)/2, sd = 5)
x <- sort_tornado(tor, bin_weights = w)
}
\seealso{
Other tornado utilities: 
\code{\link{flatten_features}()},
\code{\link{melt_tornado}()},
\code{\link{norm_tornado}()}
}
\concept{tornado utilities}
