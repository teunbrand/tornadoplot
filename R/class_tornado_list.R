# The tornado_list class is constructed with the vctrs package.
# The reason for making this class is type stability and printing niceties.
#
# An alternative I've considered is containing a CharacterList class, with
# a parallel attribute containing the dimensions. However, I ran into trouble
# shepherding these through the internals of ggplot2.
#

# Constructor -------------------------------------------------------------

#' @importFrom vctrs new_list_of new_rcrd vec_recycle
#' @importFrom rlang hash
tornado_list <- function(x = matrix(character(0), 0, 0), scale = NULL) {
  if (is.matrix(x)) {
    x <- list(x)
  }

  # Check input
  stopifnot(

    "Tornado list input should be a list." =
      is.list(x),

    "Tornado list has no matrix elements." =
      all(vapply(x, is.matrix, logical(1))),

    "Tornado list matrices are not character." =
      all(vapply(x, mode, character(1)) == "character")

  )

  if (!is.null(scale)) {
    if (inherits(scale, "ScaleContinuous")) {
      stopifnot(
        "Scale should be a valid continuous scale" =
          inherits(scale, "ScaleContinuous"),

        "Scale should have the right aesthetics" =
          "coverage" %in% scale$aesthetics
      )

      hsh <- hash(scale)
      scale <- list(scale)
      names(scale) <- hsh
    }
  } else {
    hsh <- hash(scale)
    scale <- list(NULL)
    names(scale) <- hsh
  }

  x <- new_list_of(
    x, ptype = matrix(character(), 0, 0)
  )

  if (length(hsh) %in% c(0, 1)) {
    hash <- vec_recycle(hsh, length(x))
  }

  new_tornado_list(x = x, hash = hash, scale = scale)
}

new_tornado_list <- function(
  x = new_list_of(ptype = matrix(character(), 0, 0)),
  hash = character(),
  scale = list()
) {
  new_rcrd(
    fields = list(colour_matrices = x,
                  hash = hash),
    scale_list = scale,
    class = "tornado_list"
  )
}

# Accessors ---------------------------------------------------------------

#' @importFrom vctrs field
get_tl_matrices <- function(x, unclass = FALSE) {
  if (inherits(x, "tornado_list")) {
    x <- field(x, "colour_matrices")
  }
  if (unclass) {
    unclass(x)
  } else {
    x
  }
}

get_tl_hash <- function(x) {
  field(x, "hash")
}

get_tl_scales <- function(x) {
  attr(x, "scale_list")
}

scale_list <- function(x) {
  hashes <- unique(get_tl_hash(x))
  x <- get_tl_scales(x)[hashes]
  structure(x, class = "covscale")
}

# Casting and coercion ----------------------------------------------------

# In establishing a common prototype for two tornado_lists, we need to combine
# the scale list attribute.

#' @importFrom vctrs vec_ptype2
#' @export
vec_ptype2.tornado_list.tornado_list <- function(x, y, ...) {
  scales <- c(get_tl_scales(x), get_tl_scales(y))
  new_tornado_list(scale = scales[!duplicated(names(scales))])
}

# When casting, we need the target scales but current data

#' @importFrom vctrs vec_cast
#' @export
vec_cast.tornado_list.tornado_list <- function(x, to, ...) {
  new_tornado_list(get_tl_matrices(x), get_tl_hash(x), get_tl_scales(to))
}

# The levels method is a stupid workaround for rbind

#' @export
#' @method levels tornado_list
levels.tornado_list <- function(x) {
  NULL
}

# Printing methods --------------------------------------------------------

#' @export
#' @method format tornado_list
format.tornado_list <- function(x, ...) {
  x <- get_tl_matrices(x)
  ans <- character(length(x))
  nas <- lengths(x) == 0
  ans[nas] <- NA_character_
  dim <- t(vapply(x[!nas], dim, integer(2)))
  ans[!nas] <- paste0(
    "[", format(dim[, 1]), " x ", format(dim[, 2]), "] ",
    "colour matrix"
  )
  ans
}

#' @importFrom vctrs obj_print_data
#' @export
#' @method obj_print_data tornado_list
obj_print_data.tornado_list <- function(x, ...) {
  if (length(x) == 0)
    return()
  print(format(x), quote = FALSE)
  invisible(x)
}

#' @importFrom vctrs vec_ptype_full
#' @export
#' @method vec_ptype_full tornado_list
vec_ptype_full.tornado_list <- function(x, ...) {
  "list_of_tornados"
}

#' @importFrom vctrs vec_ptype_abbr
#' @export
#' @method vec_ptype_abbr tornado_list
vec_ptype_abbr.tornado_list <- function(x, ...) {
  "lot"
}
