#' Make a tornado plot
#'
#' The `tornado_plot()` function exists to easily jot down a tornado plot. The
#' idea behind this function is that whenever you're in doubt on how to plot
#' your tornado data, this is the function that tries to figure it out for you.
#'
#' @param data A \linkS4class{TornadoExperiment}, `tornado_df` or other type of
#'   data forwarded to [`build_tornado()`][build_tornado()].
#' @param ... Arguments passed down to other functions. See the methods section
#'   below.
#' @param sort A `logical(1)`. Should the tornado be sorted from high to low?
#'   Does not apply for the `tornado_df` method.
#'
#' @return Prints a plot as a side effect. Returns a `TornadoExperiment` with
#'   the plot in the [`metadata()`][S4Vectors::Annotated-class] slot. When
#'   `data` is a `tornado_df` object, instead it returns the plot as-is.
#' @export
#' @name tornado_plot
#'
#' @examples
#' # These features and data aren't really pretty
#' feats <- dummy_features()
#' data  <- dummy_granges_data()
#'
#' # This will extract and plot the data. The extracted data is now stored in x.
#' x <- tornado_plot(data, feats, width = 3000)
#'
#' # `x` is the data underlying the plot in an TornadoExperiment object.
#' print(x)
#'
#' # The plot generated above is stored in the metadata
#' metadata(x)$plot
#'
#' # Or the plot can be recomputed from the data
#' tornado_plot(x)
setGeneric(
  "tornado_plot",
  signature = "data",
  function(data, ...) standardGeneric("tornado_plot")
)

#' @export
#' @describeIn tornado_plot The `...` argument is forwarded to
#'   [`build_tornado()`][build_tornado()].
#' @inheritDotParams build_tornado
setMethod(
  "tornado_plot",
  signature = c(data = "ANY"),
  function(data, ...) {
    data <- build_tornado(data = data, ...)
    callGeneric(data, ...)
  }
)

#' @export
#' @describeIn tornado_plot The `...` argument is forwarded to
#'   [`prep_tornado()`][prep_tornado()].
#' @inheritDotParams prep_tornado -tornado
#' @importFrom S4Vectors metadata<-
setMethod(
  "tornado_plot",
  signature = c(data = "TornadoExperiment"),
  function(data, ..., sort = NULL) {
    sort <- sort %||% metadata(data)$sorted %||% TRUE
    if (isTRUE(sort)) {
      data <- sort_tornado(data)
    }
    prep <- prep_tornado(data, ...)
    plot <- callGeneric(prep,  ...)
    metadata(data) <- list(plot = plot)
    invisible(data)
  }
)

setOldClass("tornado_df")

#' @export
#' @inheritDotParams autoplot.tornado_df -object
#' @describeIn tornado_plot The `...` argument is forwarded to
#'   [`autoplot()`][autoplot.tornado_df()].
setMethod(
  "tornado_plot",
  signature = c(data = "tornado_df"),
  function(data, ...) {
    g <- autoplot(data, ...)
    print(g)
    return(invisible(g))
  }
)
