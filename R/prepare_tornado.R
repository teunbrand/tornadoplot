#' Prepare tornado array for plotting
#'
#' Converts the data from a \linkS4class{TornadoExperiment} to a `data.frame`
#' that is accepted in ggplot2 plots.
#'
#' @param tornado A \linkS4class{TornadoExperiment} object.
#' @param assay_name A `character(1)`: one of the assay names retrieved by
#'  `assayNames(tornado)`.
#' @param upper,lower Limits to the colour scale. Can be one of the following:
#' * A `numeric` to set absolute limits directly.
#' * A `character` in the form of a number prefixed by `"p"` or `"q"` to set
#'   the percentile or quantile respectively. For example, `"p99"` and `"q0.95"`
#'   set the limit to the 99th percentile or 0.95th quantile respectively. The
#'   limits as specified in the `lower` and `upper` arguments can be overridden
#'   if the `scale` argument has non-default limits.
#' @param scale_title A `character(1)` with the title for the colourbar. Can
#'   be overridden when the `scale` argument has a non-default name.
#' @param scale A [continuous fill scale][ggplot2::scale_fill_continuous()] from
#'   the ggplot2 package to use to colour the tornado. When `NULL` (default),
#'   the scale will be retrieved from the option `"tornadoplot.default.scale"`.
#' @param ... Not currently used.
#'
#' @return A `tornado_df`, `data.frame` class object with a modified scale in
#'   the `"scale"` attribute.
#' @details The typical way of plotting a raster in ggplot2 is with the
#'   [`geom_raster()`][ggplot2::geom_raster()] function. However, this requires
#'   you to melt the raster first only to later reconstitute it, which can be
#'   rather inefficient. To circumvent this inconvenience, we pre-colour the
#'   data here with a ggplot2 scale.
#'
#' @aliases tornado_df
#' @importFrom ggplot2 guide_colourbar
#' @importFrom vctrs data_frame
#' @export
#' @md
#'
#' @examples
#' NULL
prep_tornado <- function(
  tornado, assay_name = assayNames(tornado),
  upper = "q0.99", lower = 0,
  scale_title = "Coverage",
  scale = NULL, ...
) {
  assay_name <- match.arg(assay_name, assayNames(tornado))
  data  <- assay(tornado, assay_name)

  # Resolve scale
  scale  <- resolve_scale(data, scale, upper, lower, scale_title)
  limits <- scale$get_limits()

  # Colour the data
  data <- colour_tornado(data, scale)
  # Split to list of matrices
  data <- dice_tornado(data, get_keydata(tornado, "row"))

  # Get position metadata
  nbin <- nbin(tornado)
  position <- get_bin_position(tornado)
  outerbin <- c(diff(head(position, 2)), diff(tail(position, 2)))
  position <- range(position) + c(-0.5, 0.5) * outerbin

  # Get sample metadata
  samples  <- get_sample_data(tornado)
  samples  <- as.list(samples[flat_idx(data, 2), , drop = FALSE])

  # Get feature metadata
  features <- get_keydata(tornado, "row")
  featset  <- unique(features)
  features <- runLength(Rle(features))

  # Structure
  index <- seq_along(data)

  ans <- data_frame(
    xmin = position[1],
    xmax = position[2],
    ymin = 0.5,
    ymax = features[flat_idx(data, 1)] + 0.5,
    feature_set = featset[flat_idx(data, 1)],
    .name_repair = "minimal"
  )
  cnames <- setdiff(names(samples), names(ans))
  ans[cnames] <- samples[cnames]

  # Attach heatmaps
  dim(data) <- NULL
  ans$tornado <- tornado_list(data, scale)

  # Set class and attributes
  attr(ans, "keys") <- list(sample = colKey(tornado))
  class(ans) <- c("tornado_df", class(ans))
  ans
}

#' @export
#' @importFrom vctrs vec_rbind
#' @method rbind tornado_df
rbind.tornado_df <- function(...) {
  vec_rbind(...)
}

# Helpers -----------------------------------------------------------------

## Wrangling --------------------------------------------------------------

split_array_rows <- function(x, f, drop = FALSE, ...) {
  if (length(dim(x)) == 2) {
    fun <- function(i) x[i, , drop = FALSE]
  } else {
    fun <- function(i) x[i, , , drop = FALSE]
  }
  lapply(split(x = seq_len(nrow(x)), f, drop = drop, ...), fun)
}

# Convert array of number to array of colours
colour_tornado <- function(tornado, scale) {
  # Store and remove dimensions
  old_dim <- dim(tornado)
  dim(tornado) <- NULL

  # Get limits
  limits <- scale$get_limits()

  # Squish out-of-bounds values
  tornado <- scale$trans$transform(tornado)
  tornado <- pmin(pmax(tornado, limits[1]), limits[2])
  tornado <- scale$rescale(tornado)

  # Apply palette
  uniq     <- unique(tornado)
  pal      <- scale$palette(uniq)
  tornado  <- pal[match(tornado, uniq)]

  # Reset dimensions
  dim(tornado) <- old_dim
  tornado
}

## Scales ------------------------------------------------------------------

#' @importFrom ggplot2 waiver guide_colourbar
resolve_scale <- function(
  data,
  scale = NULL,
  upper,
  lower,
  title = NULL
) {
  if (is.null(scale)) {
    scale <- getOption(
      "tornadoplot.default.scale",
      default = tornado_default_scale()
    )
  }
  if (!inherits(scale, "ScaleContinuous")) {
    stop("Cannot interpret scale as continuous colour scale.")
  }

  # Set scale attributes
  scale$aesthetics <- "coverage"

  # Set limits
  limits <- choose_limits(
    lower  = lower, upper = upper,
    oldlim = scale$limits,
    data   = data,
    trans  = scale$trans$transform
  )
  scale$limits <- limits

  # Set colourbar guide
  name <- scale$name
  if (inherits(name, "waiver")) {
    name <- title
  }
  guide <- scale$guide
  if (is.character(guide)) {
    # Should fail here if there is no such function
    guide <- match.fun(paste0("guide_", guide))
  }
  if (is.function(guide)) {
    guide <- guide()
  }
  if (inherits(guide, "guide")) {
    guide$available_aes <- "coverage"
    guide$title <- title
  }
  scale$guide <- guide

  # Finally, train data
  scale$train(data)
  scale
}

# Parses text limits when appropriate
choose_limits <- function(lower = 0, upper = "q0.99",
                          oldlim = NULL, data = c(0, 1),
                          trans = force) {
  if (is.null(oldlim)) {
    # Determine upper bound of colour scale
    if (!is.numeric(upper)) {
      if (grepl("^p|^q", upper)) {
        num_upper <- as.numeric(gsub("^p|^q", "", upper))
        if (grepl("^p", upper)) {
          num_upper <- num_upper / 100
        }
        upper <- quantile(data, num_upper)
      } else {
        warning("Could not interpret 'upper'.")
        upper <- NA
      }
    }
    # Determine lower bound of colour scale
    if (!is.numeric(lower)) {
      if (grepl("^p|^q", lower)) {
        num_lower <- as.numeric(gsub("^p|^q", "", lower))
        if (grepl("^p", lower)) {
          num_lower <- num_lower / 100
        }
        lower <- quantile(data, num_lower)
      } else {
        warning("Could not interpret 'lower'.")
        lower <- NA
      }
    }
    limits <- trans(unname(sort(c(lower, upper))))
  } else {
    limits <- oldlim
  }
  return(limits)
}

# Colours are from the following code:
# colorspace::sequential_hcl(12, palette = "PuBu", rev = TRUE)
# But the colours are hardcoded to avoid a dependency
tornado_default_scale <- function(...) {
  ggplot2::scale_fill_gradientn(
    colours = c("#F8F9FF", "#EBEDF7", "#DADEED", "#C6CDE4",
                "#B0BCDB", "#97ABD2", "#7A99CA", "#5788C1",
                "#1577B9", "#06659A", "#06527B", "#0E3F5C"),
    name = "Coverage",
    ...
  )
}
