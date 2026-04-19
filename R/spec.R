# grass_spec S3 class and family-specific constructors.
#
# The spec tells grass_report() which metric family to compute and which
# reference curve to attach. Today only grass_spec_binary() returns a
# usable spec; the ordinal / multirater / continuous constructors return
# placeholder specs that error with a pointer to ?grass_roadmap when
# they reach a dispatch generic.

new_grass_spec <- function(family, ...) {
  params <- list(...)
  structure(
    c(list(family = family), params),
    class = c(paste0("grass_spec_", family), "grass_spec")
  )
}

# ---- Binary (implemented) ---------------------------------------------

#' Binary two-rater agreement spec
#'
#' The spec for the currently-implemented metric family: Cohen's kappa,
#' PABAK, and Gwet's AC1 for two raters on a binary outcome, evaluated
#' against a prevalence-conditioned reference curve at the chosen
#' `reference_level`.
#'
#' @param reference_level One of `0.70`, `0.80`, `0.85`, `0.90`. Selects the
#'   Se = Sp rater-quality band whose analytical Youden-J reference curve
#'   is attached to the result. `0.85` is the default and reproduces the
#'   "high" quality curve from prior releases; `0.70` reproduces the
#'   "medium" curve. `0.80` and `0.90` are new in 0.1.2. Pass `NULL` to
#'   skip the reference attachment entirely (metrics in context only).
#'
#' @return A `grass_spec_binary` object carrying `$family = "binary"` and
#'   `$reference_level`.
#' @export
#'
#' @examples
#' grass_spec_binary()
#' grass_spec_binary(reference_level = 0.80)
grass_spec_binary <- function(reference_level = 0.85) {
  valid <- c(0.70, 0.80, 0.85, 0.90)
  if (!is.null(reference_level)) {
    if (!is.numeric(reference_level) || length(reference_level) != 1 ||
        !any(abs(reference_level - valid) < 1e-8)) {
      stop("`reference_level` must be one of 0.70, 0.80, 0.85, 0.90, or NULL. ",
           "Got: ", deparse(reference_level), ".", call. = FALSE)
    }
    reference_level <- valid[which.min(abs(reference_level - valid))]
  }
  new_grass_spec("binary", reference_level = reference_level)
}

# ---- Stub constructors (error at dispatch) ----------------------------

#' Ordinal two-rater agreement spec (planned)
#'
#' Placeholder for the GRASS ecosystem's ordinal family: weighted Cohen's
#' kappa (linear and quadratic weights), Gwet's AC2, and ordinal
#' Krippendorff's alpha, each with prevalence-conditioned reference
#' curves analogous to the binary family.
#'
#' Constructing the spec is legal (so users can write code that is ready
#' for the release); passing it to [grass_report()] errors until the
#' implementation lands. See `?grass_roadmap` for the framework taxonomy.
#'
#' @param ... Ignored in 0.1.2; reserved for future parameters.
#' @return A placeholder `grass_spec_ordinal` object.
#' @seealso [grass_roadmap]
#' @export
grass_spec_ordinal <- function(...) {
  new_grass_spec("ordinal")
}

#' Multi-rater nominal agreement spec (planned)
#'
#' Placeholder for Fleiss' kappa and nominal Krippendorff's alpha (more
#' than two raters, nominal scale), with prevalence-conditioned
#' reference curves analogous to the binary family.
#'
#' See `?grass_roadmap`.
#'
#' @param ... Ignored in 0.1.2; reserved for future parameters.
#' @return A placeholder `grass_spec_multirater` object.
#' @seealso [grass_roadmap]
#' @export
grass_spec_multirater <- function(...) {
  new_grass_spec("multirater")
}

#' Continuous agreement spec (planned)
#'
#' Placeholder for the continuous family: intraclass correlation
#' coefficient (Shrout-Fleiss forms), Lin's concordance correlation
#' coefficient, Bland-Altman limits of agreement (parametric and
#' non-parametric), standard error of measurement, and generalizability
#' theory variance components.
#'
#' See `?grass_roadmap`.
#'
#' @param ... Ignored in 0.1.2; reserved for future parameters.
#' @return A placeholder `grass_spec_continuous` object.
#' @seealso [grass_roadmap]
#' @export
grass_spec_continuous <- function(...) {
  new_grass_spec("continuous")
}

# ---- Predicate and print helpers --------------------------------------

is_grass_spec <- function(x) inherits(x, "grass_spec")

#' @export
print.grass_spec <- function(x, ...) {
  cat("<grass_spec> family = ", x$family, "\n", sep = "")
  params <- x[setdiff(names(x), "family")]
  if (length(params)) {
    for (nm in names(params)) {
      val <- params[[nm]]
      cat("  ", nm, " = ",
          if (is.null(val)) "NULL" else format(val), "\n", sep = "")
    }
  }
  invisible(x)
}
