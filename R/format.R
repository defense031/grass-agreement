#' One-line paper-ready summary of a grass report
#'
#' Returns a single character string suitable for dropping into a
#' manuscript line or a results table cell. Contains the three agreement
#' metrics (with a confidence interval on kappa), the skew diagnostics
#' (PI, BI), sample size, and regime.
#'
#' @param x A `grass_result` object.
#' @param digits Number of decimal places for rounding. Default 2.
#' @param ascii If `TRUE` (the default), emit `kappa` instead of the Unicode
#'   `\u03ba`. Safer for Slack, Markdown, and non-UTF-8 locales. Pass
#'   `ascii = FALSE` for a Unicode-rendering context (e.g., RStudio console).
#' @param ci_width If `TRUE`, append the kappa CI width and a one-word
#'   descriptor (`tight` / `moderate` / `wide`) based on Wilson-logit
#'   half-width cutpoints at 0.10 and 0.20. These cutpoints are calibrated
#'   for kappa CIs; Se/Sp CIs at the same N typically run roughly 2-2.5x
#'   tighter. Default `FALSE`.
#'
#' @return A single character string.
#' @export
#'
#' @examples
#' tab <- matrix(c(88, 10, 14, 88), nrow = 2,
#'               dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
#' grass_format_report(grass_report(tab, format = "matrix"))
#' grass_format_report(grass_report(tab, format = "matrix"), ci_width = TRUE)
grass_format_report <- function(x, digits = 2, ascii = TRUE, ci_width = FALSE) {
  if (!inherits(x, "grass_result")) {
    stop("`x` must be a grass_result, as returned by grass_report().",
         call. = FALSE)
  }
  v <- x$metrics$values
  kappa_name <- if (isTRUE(ascii)) "kappa" else "\u03ba"
  fmt <- paste0("%.", digits, "f")
  ci  <- sprintf(paste0("[", fmt, ", ", fmt, "]"),
                 unname(v["kappa_wilson_lower"]),
                 unname(v["kappa_wilson_upper"]))
  base <- sprintf(
    paste0("%s = ", fmt, " ", "%s, PABAK = ", fmt, ", AC1 = ", fmt,
           ", PI = ", fmt, ", BI = ", fmt, ", N = %d, prevalence = ", fmt,
           ", %s regime"),
    kappa_name,
    unname(v["kappa"]),
    ci,
    unname(v["PABAK"]),
    unname(v["AC1"]),
    unname(v["prevalence_index"]),
    unname(v["bias_index"]),
    x$metrics$n,
    x$prevalence,
    x$regime
  )
  if (isTRUE(ci_width)) {
    lo <- unname(v["kappa_wilson_lower"])
    hi <- unname(v["kappa_wilson_upper"])
    width <- hi - lo
    # Round to display precision before comparing. (0.50 - 0.40)/2 computes
    # to 0.0499999... in binary float; a naive `half < 0.10` flips the
    # descriptor at exact cutpoints. Rounding to 4 decimals matches the
    # precision the user sees in the printed interval and removes the
    # boundary-flip surprise.
    half  <- round(width / 2, 4)
    # Cutpoints calibrated for kappa Wilson-logit CIs, which run roughly
    # 2-2.5x wider than component Se/Sp CIs at the same N. The 0.10 / 0.20
    # half-width thresholds match empirical kappa CI spread rather than the
    # Se/Sp conventions (0.05 / 0.10) from clinical diagnostic accuracy.
    descriptor <- if (is.na(half))       "unknown"
                  else if (half <  0.10) "tight"
                  else if (half <  0.20) "moderate"
                  else                   "wide"
    base <- paste0(base,
                   sprintf(paste0("; kappa CI width = ", fmt, ", %s"),
                           width, descriptor))
  }
  base
}
