# Internal helpers.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Emit a message at most once per R session, keyed by a short string.
# Uses the package environment initialized in zzz.R.
msg_once <- function(key, msg) {
  if (is.null(.grass_env$msg_seen[[key]])) {
    message(msg)
    .grass_env$msg_seen[[key]] <- TRUE
  }
  invisible()
}

#' Reset grass's once-per-session message cache
#'
#' `grass` emits coercion messages and keyword-fallthrough warnings at most
#' once per R session, keyed by the column name. When looping
#' `grass_report()` over many studies with consistently named columns, only
#' the first study in the loop sees the message. Call `reset_grass_warnings()`
#' to clear the cache so the next call re-evaluates.
#'
#' @return `invisible(NULL)`.
#' @export
#'
#' @examples
#' reset_grass_warnings()
reset_grass_warnings <- function() {
  .grass_env$msg_seen <- list()
  invisible(NULL)
}

# Validate prevalence is a single finite numeric in (0, 1). Clamp values in
# (0, 0.01) and (0.99, 1) with a warning; exact 0 and 1 are degenerate (no
# agreement problem exists to compare against a reference) and error.
validate_prevalence <- function(p) {
  if (!is.numeric(p) || length(p) != 1 || !is.finite(p)) {
    stop("`prevalence` must be a single finite numeric value.", call. = FALSE)
  }
  if (p < 0 || p > 1) {
    stop("`prevalence` must lie in [0, 1]. Got ", p, ".", call. = FALSE)
  }
  if (p == 0 || p == 1) {
    stop("`prevalence` = ", p,
         " is degenerate (all subjects in one class). No reference comparison ",
         "is defined. Supply a value in (0, 1) or estimate from data.",
         call. = FALSE)
  }
  if (p < 0.01) {
    warning("Prevalence ", p, " is below 0.01; clamping to 0.01 for reference lookup.",
            call. = FALSE)
    return(0.01)
  }
  if (p > 0.99) {
    warning("Prevalence ", p, " is above 0.99; clamping to 0.99 for reference lookup.",
            call. = FALSE)
    return(0.99)
  }
  p
}

# Pretty-print a numeric (scalar or vector) to 4 decimal places; NA -> "NA".
fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "   NA", formatC(x, format = "f", digits = digits))
}

# Pretty-print a confidence interval.
fmt_ci <- function(lower, upper, digits = 4) {
  if (is.na(lower) || is.na(upper)) return("[NA, NA]")
  paste0("[", fmt_num(lower, digits), ", ", fmt_num(upper, digits), "]")
}
