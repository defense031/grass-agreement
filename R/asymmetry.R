# Column A of the GRASS Reporting Card: per-rater Se/Sp asymmetry diagnostic.
# Assigns a model-safety tier (ok / caution / unsafe) from a scalar δ̂
# summary of per-rater |Se − Sp| gaps. Ships with grass 0.2.x; see
# paper2/review/framework_notes.md §0.6 for the three-tier architecture
# this function operationalises.

#' Column A of the Reporting Card: model-safety tier from rater asymmetry
#'
#' `check_asymmetry()` takes per-rater sensitivity and specificity estimates,
#' computes the per-rater gap `|Ŝe_j − Ŝp_j|`, reduces them to a scalar
#' `δ̂`, and assigns the model-safety tier that governs whether `q̂` (the
#' GRASS operating-quality projection onto the Se = Sp diagonal) is
#' trustworthy as the primary summary. This is the Column A output of the
#' §8 Reporting Card in the merged GRASS binary-rater-reliability paper.
#'
#' Tier thresholds follow the three-tier architecture:
#'
#' - **Tier 1 — `ok`** (`δ̂ < 0.05`): diagonal default. Report `q̂ ± SE`
#'   plus EDR and EMR_panel without per-rater disclosure.
#' - **Tier 2 — `caution`** (`0.05 ≤ δ̂ < 0.10`): diagonal + diagnostic.
#'   Report `q̂ ± SE` with a caution flag and the `δ̂` value; EDR / EMR_panel
#'   carry a footnote that rater asymmetry is non-negligible.
#' - **Tier 3 — `unsafe`** (`δ̂ ≥ 0.10`): full latent-class. `q̂` is withheld
#'   as primary; report per-rater `(Ŝe_j, Ŝp_j)` via a Hui-Walter / Dawid-
#'   Skene fit. PABAK's prevalence-flatness was conditional on Se = Sp; at
#'   Tier 3 that condition is visibly broken.
#'
#' `se` and `sp` are not identifiable from a single 2-rater binary table
#' without external ground truth. Supply them from: (a) a simulation with
#' known truth, (b) a Hui-Walter / Dawid-Skene latent-class fit under
#' `k ≥ 3` and prevalence heterogeneity, or (c) a reference-standard
#' comparison. Passing raw `table()` off-diagonals as `se` / `sp` is
#' **wrong** and the function cannot detect the error.
#'
#' @param se Numeric vector, one per rater. Per-rater sensitivity estimates
#'   in `[0, 1]`.
#' @param sp Numeric vector, one per rater. Per-rater specificity estimates
#'   in `[0, 1]`. Same length as `se`.
#' @param rater Optional character vector of rater labels, same length as
#'   `se`. Defaults to `"R1"`, `"R2"`, ... .
#' @param threshold_caution Boundary between Tier 1 (`ok`) and Tier 2
#'   (`caution`). Default `0.05`.
#' @param threshold_unsafe Boundary between Tier 2 (`caution`) and Tier 3
#'   (`unsafe`). Default `0.10`.
#' @param summary How to reduce per-rater gaps to the scalar `δ̂`. `"max"`
#'   (default, conservative tripwire — any single rater above the threshold
#'   triggers escalation) or `"mean"` (panel-average asymmetry). The
#'   framework uses `"max"`; `"mean"` is provided for sensitivity checks.
#'
#' @return An S3 object of class `grass_asymmetry` with fields:
#' - `per_rater`: data.frame of `rater`, `se`, `sp`, and `gap = |se − sp|`
#' - `delta_hat`: the scalar `δ̂` summary
#' - `summary`: which summary statistic was used (`"max"` or `"mean"`)
#' - `tier`: integer tier (`1`, `2`, or `3`)
#' - `regime`: character regime label (`"ok"`, `"caution"`, or `"unsafe"`)
#' - `thresholds`: named list of the caution and unsafe cutoffs
#'
#' @seealso [classify()] for the companion Column B tier (use-case fit on
#'   EMR_panel).
#' @export
#'
#' @examples
#' # Tier 1: symmetric raters — safe to use q̂ as primary
#' check_asymmetry(se = c(0.86, 0.88, 0.84), sp = c(0.85, 0.87, 0.86))
#'
#' # Tier 2: one rater pressing the Se = Sp assumption
#' check_asymmetry(se = c(0.90, 0.88, 0.92), sp = c(0.82, 0.86, 0.85))
#'
#' # Tier 3: within-rater Se-favoring regime — requires latent-class fit
#' check_asymmetry(se = c(0.95, 0.93, 0.94), sp = c(0.78, 0.80, 0.79))
check_asymmetry <- function(se, sp,
                            rater = NULL,
                            threshold_caution = 0.05,
                            threshold_unsafe = 0.10,
                            summary = c("max", "mean")) {
  summary <- match.arg(summary)

  if (!is.numeric(se) || !is.numeric(sp)) {
    stop("`se` and `sp` must be numeric vectors.", call. = FALSE)
  }
  if (length(se) != length(sp)) {
    stop("`se` and `sp` must have the same length (one entry per rater).",
         call. = FALSE)
  }
  if (length(se) < 1L) {
    stop("Need at least one rater's (se, sp) estimate.", call. = FALSE)
  }
  if (any(!is.finite(se)) || any(!is.finite(sp))) {
    stop("`se` and `sp` must be finite.", call. = FALSE)
  }
  if (any(se < 0 | se > 1) || any(sp < 0 | sp > 1)) {
    stop("`se` and `sp` values must lie in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(threshold_caution) || !is.numeric(threshold_unsafe) ||
      length(threshold_caution) != 1L || length(threshold_unsafe) != 1L ||
      threshold_caution <= 0 || threshold_unsafe <= threshold_caution) {
    stop("Thresholds must satisfy 0 < threshold_caution < threshold_unsafe.",
         call. = FALSE)
  }

  if (is.null(rater)) {
    rater <- paste0("R", seq_along(se))
  } else if (length(rater) != length(se)) {
    stop("`rater` must be the same length as `se` and `sp`.", call. = FALSE)
  }

  gap <- abs(se - sp)
  delta_hat <- if (summary == "max") max(gap) else mean(gap)

  tier <- if (delta_hat < threshold_caution) {
    1L
  } else if (delta_hat < threshold_unsafe) {
    2L
  } else {
    3L
  }
  regime <- c("ok", "caution", "unsafe")[tier]

  out <- list(
    per_rater = data.frame(
      rater = as.character(rater),
      se = as.numeric(se),
      sp = as.numeric(sp),
      gap = as.numeric(gap),
      stringsAsFactors = FALSE
    ),
    delta_hat = as.numeric(delta_hat),
    summary = summary,
    tier = tier,
    regime = regime,
    thresholds = list(caution = as.numeric(threshold_caution),
                      unsafe = as.numeric(threshold_unsafe))
  )
  class(out) <- c("grass_asymmetry", "list")
  out
}

#' @export
print.grass_asymmetry <- function(x, digits = 3, ...) {
  cat("grass asymmetry diagnostic (Column A of Reporting Card)\n",
      "  raters              : ", nrow(x$per_rater), "\n",
      "  per-rater gaps |Se-Sp|\n", sep = "")
  pr <- x$per_rater
  for (i in seq_len(nrow(pr))) {
    cat(sprintf("    %-8s  Se = %.*f   Sp = %.*f   gap = %.*f\n",
                pr$rater[i],
                digits, pr$se[i],
                digits, pr$sp[i],
                digits, pr$gap[i]))
  }
  cat(sprintf("  delta_hat (%s)    : %.*f\n",
              x$summary, digits, x$delta_hat))
  cat(sprintf("  thresholds           : caution = %.*f,  unsafe = %.*f\n",
              digits, x$thresholds$caution,
              digits, x$thresholds$unsafe))
  cat("  tier                 : ", x$tier,
      "  (", x$regime, ")\n", sep = "")
  tier_msg <- switch(
    x$regime,
    ok      = "Diagonal default: report q_hat plus EDR / EMR_panel.",
    caution = "Diagonal + diagnostic: report q_hat with caution flag and delta_hat.",
    unsafe  = "Full latent-class: withhold q_hat as primary; report per-rater (Se_j, Sp_j)."
  )
  cat("  reporting guidance   : ", tier_msg, "\n", sep = "")
  invisible(x)
}

#' Coerce a `grass_asymmetry` result to a one-row data.frame
#'
#' @param x A `grass_asymmetry` object.
#' @param row.names,optional,... Standard arguments; ignored except
#'   `row.names`.
#' @return A one-row data.frame with `n_raters`, `delta_hat`, `summary`,
#'   `tier`, and `regime` columns.
#' @export
as.data.frame.grass_asymmetry <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    n_raters = nrow(x$per_rater),
    delta_hat = x$delta_hat,
    summary = x$summary,
    tier = x$tier,
    regime = x$regime,
    threshold_caution = x$thresholds$caution,
    threshold_unsafe = x$thresholds$unsafe,
    row.names = row.names,
    stringsAsFactors = FALSE
  )
}
