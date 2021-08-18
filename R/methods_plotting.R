# Tornadoplot -------------------------------------------------------------

#' Make a tornado plot.
#'
#' The `tornado_plot()` function exists to easily jot down a tornado plot. It
#' will print the plot and return the data used for plotting with the plot
#' itself in the metadata. This function performs the tasks of retrieving and
#' aggregating data from its source with [`build_tornado()`], optionally sorting
#' the result with [`sort_tornado`] and finally plotting.
#'
#' @inheritDotParams build_tornado
#' @param sort A `logical(1)`. Should the tornado be sorted from high to low?
#'
#' @return A \linkS4class{TornadoExperiment} object with a `ggplot` object in
#'   the [`metadata()`][S4Vectors::Annotated-class] annotation.
#' @export
#' @importFrom S4Vectors metadata<-
#'
#' @examples
#' feats <- dummy_features()
#' file  <- dummy_granges_data()
#'
#' tornado_plot(feats, file, width = 3000)
tornado_plot <- function(..., sort = TRUE) {
  tor <- build_tornado(...)
  if (isTRUE(sort)) {
    tor <- sort_tornado(tor)
  }
  df <- prep_tornado(tor)
  g  <- autoplot(df)
  print(g)
  metadata(tor) <- list(plot = g)
  invisible(tor)
}

# Autolayer ---------------------------------------------------------------

#' @importFrom ggplot2 autolayer
#' @method autolayer tornado_df
#' @export
autolayer.tornado_df <- function(object, ...) {
  vars <- c("x", "y")
  vars <- c(paste0(vars, "min"), paste0(vars, "max"),
            "tornado")
  if (!all(vars %in% colnames(object))) {
    stop("Cannot find all relevant variables. Has the tornado_df been
         altered?")
  }

  list(
    geom_tornado(
      data = object,
      aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax,
        tornado  = .data$tornado,
      ),
      ...
    ),
    get_tl_scales(object$tornado)[[1]]
  )
}

#' @importFrom ggplot2 autoplot ggplot
#' @method autoplot tornado_df
#' @export
autoplot.tornado_df <- function(
  object, ...,
  facet = TRUE, x_scale = TRUE, y_scale = TRUE
) {
  g <- ggplot() +
    autolayer(object, ...)

  # Set axis titles manually instead of via scales.
  # This prevents them being 'baked in',
  # so they can be overridden with `+ labs(x = ..., y = ...)`.
  g$labels[["x"]] <- "Position relative to center"
  g$labels[["y"]] <- "Features"

  if (isTRUE(x_scale)) {
    x_scale <-
      ggplot2::scale_x_continuous(
        expand = c(0, 1),
        breaks = function(x) {
          dodge  <- diff(range(x)) * 0.1
          breaks <- scales::extended_breaks()(x)
          breaks[breaks > (x[1] + dodge) & breaks < (x[2] - dodge)]
        }
      )
  }
  if (inherits(x_scale, "ScaleContinuousPosition")) {
    g <- g + x_scale
  }
  if (isTRUE(y_scale)) {
    y_lim <- range(object$ymin, object$ymax, object$y)
    y_breaks <- scales::extended_breaks()(y_lim)
    y_breaks <- sort(union(y_breaks, y_lim - c(0, 0.5)))
    y_scale <-
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(0.5, NA),
        breaks = y_breaks,
        labels = function(x) {
          c(rep("", length(x) - 1), x[length(x)])
        }
      )
  }
  if (inherits(y_scale, "ScaleContinuousPosition")) {
    g <- g + y_scale
  }
  if (isTRUE(facet)) {
    g <- g + facet_tornado(object)
  }
  g
}

#' @method autoplot TornadoExperiment
#' @export
autoplot.TornadoExperiment <- function(
  object, ...,
  upper = "q0.99", lower = 0,
  scale_title = "Coverage", scale = NULL
) {
  object <- prep_tornado(object, upper = upper, lower = lower,
                         scale_title = scale_title, scale = scale)
  autoplot(object, ...)
}


# Geoms -------------------------------------------------------------------

## Constructor ------------------------------------------------------------

#' Tornado Heatmap
#'
#' `geom_tornado()` serves as a plotting layer for tornado heatmaps.
#'
#' @inheritParams ggplot2::geom_raster
#' @param data A `tornado_df` data.frame, such as one that results from
#'   processing a tornado array with [`prep_tornado()`][prep_tornado()].
#'
#' @details The required aesthetics are synonymous with the column names of a
#'  default `tornado_df`.
#'
#' @section Aesthetics:
#' `geom_tornado` understands the following aesthetics (required aesthetics
#' are in bold):
#'
#' * \strong{xmin}
#' * \strong{xmax}
#' * \strong{ymin}
#' * \strong{ymax}
#' * \strong{tornado}
#'
#' @return A `GeomTornado` ggproto object that can be added to a plot.
#' @export
#' @md
#' @importFrom ggplot2 layer
#'
#' @examples
#' NULL
geom_tornado <- function(
  mapping = NULL,
  data = NULL,
  position = "identity",
  ...,
  interpolate = FALSE,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  # We use a dummy for coverage so that the guide_geom method picks up on
  # shared aesthetics
  if (!("coverage" %in% names(mapping))) {
    mapping[["coverage"]] <- NA_real_
  }
  layer(
    data = data,
    mapping = mapping,
    stat = "identity",
    geom = GeomTornado,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      interpolate = interpolate,
      na.rm = na.rm,
      ...
    )
  )
}

## ggproto class ----------------------------------------------------------

#' @importFrom ggplot2 ggproto draw_key_rect aes Geom
#' @importFrom grid rasterGrob grobTree rectGrob gpar
GeomTornado <- ggproto(
  "GeomTornado",
  Geom,
  required_aes = c("xmin", "xmax", "ymin", "ymax", "tornado"),
  optional_aes = c("coverage"),
  default_aes  = aes(coverage = "grey20"),
  draw_panel = function(data, panel_params, coord, interpolate = FALSE) {

    if (!inherits(coord, "CoordCartesian")) {
      stop("geom_tornado only works with cartesian coordinates.")
    }

    data <- coord$transform(data, panel_params)

    x <- (data$xmin + data$xmax) / 2
    y <- (data$ymin + data$ymax) / 2
    width  <- abs(data$xmax - data$xmin)
    height <- abs(data$ymax - data$ymin)

    grobs <- mapply(
      function(x, y, width, height, tornado) {
        rasterGrob(
          image = tornado,
          x = x, y = y, width = width, height = height,
          default.units = "native", interpolate = interpolate
        )
      },
      x = x, y = y, width = width, height = height,
      tornado  = get_tl_matrices(data$tornado, TRUE),
      SIMPLIFY = FALSE
    )
    grobs <- do.call(grobTree, grobs)
    return(grobs)
  },
  draw_key = function(data, params, size) {
    if (is.atomic(data)) {
      data <- data.frame(fill = data)
    }
    colnames(data)[colnames(data) == "coverage"] <- "fill"
    data$alpha <- NA
    draw_key_rect(data, params, size)
  }
)

# Accessories -------------------------------------------------------------

#' Automatic facets for tornado heatmaps
#'
#' @param tornado A `tornado_df` data.frame, such as one that results from
#'   processing a tornado array with [`prep_tornado()`][prep_tornado()].
#' @inheritParams ggplot2::facet_grid
#' @inheritParams ggplot2::theme
#'
#' @return A `FacetGrid` ggproto object that can be added to a plot.
#' @export
#' @importFrom rlang .data
#' @importFrom ggplot2 facet_grid vars theme facet_null
#'
#' @examples
#' NULL
facet_tornado <- function(
  tornado,
  rows = NULL,
  cols = NULL,
  scales = NULL,
  space = NULL,
  switch = NULL,
  shrink = NULL,
  labeller = NULL,
  as.table = NULL,
  drop = NULL,
  margins = NULL,
  strip.placement = NULL
) {

  if (is.null(cols)) {
    samp <- attr(tornado, "keys")$samp %||% "sample_name"
    if (!(samp %in% names(tornado))) {
      samp <- colnames(tornado)[5]
      message("Don't know what should be considered a sample.")
    }
    sampdata <- tornado[[samp]]
    if (length(unique(sampdata)) > 1) {
      cols <- vars(.data[[samp]])
    } else {
      cols <- NULL
    }
  }
  if (is.null(rows)) {
    if (length(unique(tornado$feature_set)) > 1) {
      rows <- vars(.data$feature_set)
    } else {
      rows <- NULL
    }
  }
  theme <- theme(strip.placement = strip.placement %||% "outside")
  if (!is.null(rows) || !is.null(cols)) {
    if (!is.null(cols)) {
      facet <- facet_grid(
        rows = rows, cols = cols,
        scales = scales %||% "free_y",
        space  = space  %||% "free_y",
        shrink = shrink %||% TRUE,
        switch = switch %||% "y",
        labeller = labeller %||% "label_value",
        as.table = as.table %||% TRUE,
        drop = drop %||% TRUE,
        margins = margins %||% FALSE
      )
    } else {
      facet <- facet_grid(
        rows = rows, cols = cols,
        scales = scales %||% "fixed",
        space  = space  %||% "fixed",
        shrink = shrink %||% TRUE,
        switch = switch,
        labeller = labeller %||% "label_value",
        as.table = as.table %||% TRUE,
        drop = drop %||% TRUE,
        margins = margins %||% FALSE
      )
    }
    facet <- list(facet, theme)
  } else {
    facet <- facet_null(shrink = shrink %||% TRUE)
  }
  return(facet)
}

#' @export
#' @importFrom ggplot2 scale_type
#' @method scale_type tornado_list
scale_type.tornado_list <- function(x) {
  "identity"
}


