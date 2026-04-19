#' @export
summary.grass_metrics <- function(object, ...) {
  structure(object, class = c("summary.grass_metrics", class(object)))
}

#' @export
print.summary.grass_metrics <- function(x, digits = 4, ...) {
  v <- x$values
  cat("grass metrics -- summary\n")
  cat("  N = ", x$n, "   positive level = ", shQuote(x$positive_level), "\n", sep = "")
  cat("  2x2 table\n")
  print(x$table)
  cat("\n")
  cat("  Agreement\n")
  cat("    P0 (observed)      : ", fmt_num(v["P0"], digits), "\n", sep = "")
  cat("    Pe (expected)      : ", fmt_num(v["Pe"], digits), "\n", sep = "")
  cat("    Cohen's kappa      : ", fmt_num(v["kappa"], digits), "\n", sep = "")
  cat("      Wald 95% CI      : ",
      fmt_ci(v["kappa_wald_lower"], v["kappa_wald_upper"], digits), "\n", sep = "")
  cat("      Wilson-logit 95% : ",
      fmt_ci(v["kappa_wilson_lower"], v["kappa_wilson_upper"], digits), "\n", sep = "")
  cat("    PABAK              : ", fmt_num(v["PABAK"], digits), "\n", sep = "")
  cat("    Gwet's AC1         : ", fmt_num(v["AC1"], digits), "\n", sep = "")
  cat("    Positive agreement : ", fmt_num(v["pos_agreement"], digits), "\n", sep = "")
  cat("    Negative agreement : ", fmt_num(v["neg_agreement"], digits), "\n", sep = "")
  cat("\n  Skew diagnostics\n")
  cat("    Prevalence index   : ", fmt_num(v["prevalence_index"], digits), "\n", sep = "")
  cat("    Bias index         : ", fmt_num(v["bias_index"], digits), "\n", sep = "")
  invisible(x)
}

#' @export
summary.grass_result <- function(object, ...) {
  structure(object, class = c("summary.grass_result", class(object)))
}

#' @export
print.summary.grass_result <- function(x, digits = 4, ...) {
  print.grass_result(x, digits = digits)
  invisible(x)
}

#' @export
summary.grass_reference <- function(object, ...) {
  structure(object, class = c("summary.grass_reference", class(object)))
}

#' @export
print.summary.grass_reference <- function(x, digits = 4, ...) {
  print.grass_reference(x, digits = digits)
  invisible(x)
}
