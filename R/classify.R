# Column B of the GRASS Reporting Card: use-case tier on EMR_panel.
# Assigns one of Fails / Indeterminate / Meets / Met with distinction /
# Exceeds against a declared tolerance T on the 95% CI of EMR_panel. Ships
# with grass 0.2.x. See paper2/review/framework_notes.md §0.6.5 for the
# label system and why it does not re-create Landis-Koch.

# ---- Use-case ladder --------------------------------------------------

# Illustrative defaults. The MEADOW Field Guide carries the canonical
# table; practitioners may declare custom tolerances via `tolerance = `.
.grass_use_case_ladder <- list(
  research   = list(T = 0.15, source = "Research coding with adjudication (low-stakes)"),
  screening  = list(T = 0.10, source = "Clinical screening (USPSTF-style)"),
  diagnostic = list(T = 0.05, source = "Clinical test validation (CLSI-style)"),
  clinical   = list(T = 0.02, source = "High-stakes clinical decision (FDA qualified-biomarker)")
)

# Rank order (most permissive first); used to locate T_next (one step
# harder on the ladder).
.grass_use_case_order <- c("research", "screening", "diagnostic", "clinical")

#' Return the GRASS use-case tolerance ladder
#'
#' Returns the illustrative defaults from the merged GRASS binary-rater-
#' reliability paper. The MEADOW Field Guide is the canonical source;
#' practitioners may declare custom tolerances per study with explicit
#' disclosure.
#'
#' @return A data.frame with columns `use_case`, `tolerance`, and `source`.
#' @export
#'
#' @examples
#' grass_use_case_ladder()
grass_use_case_ladder <- function() {
  data.frame(
    use_case = .grass_use_case_order,
    tolerance = vapply(.grass_use_case_order,
                       function(nm) .grass_use_case_ladder[[nm]]$T,
                       numeric(1)),
    source = vapply(.grass_use_case_order,
                    function(nm) .grass_use_case_ladder[[nm]]$source,
                    character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# Locate T_next for a named use case. Returns NA_real_ at the top of the
# ladder (no harder use case exists).
.lookup_tolerance_next <- function(use_case) {
  idx <- match(use_case, .grass_use_case_order)
  if (is.na(idx) || idx >= length(.grass_use_case_order)) {
    return(NA_real_)
  }
  nxt <- .grass_use_case_order[idx + 1L]
  .grass_use_case_ladder[[nxt]]$T
}

# ---- EMR_panel helper --------------------------------------------------

#' Panel-level misclassification rate under k-rater majority vote
#'
#' Under the Se = Sp = q diagonal (GRASS Tier 1), the panel's majority
#' vote across `k` raters agrees with the latent class with probability
#' `P(Bin(k, q) >= ceiling(k/2))`. `emr_panel()` returns `1 -` that
#' probability, which is the rate at which the panel's decision
#' misclassifies a randomly drawn subject.
#'
#' Only odd `k` is supported in this release; even `k` requires a
#' tie-breaking rule that the framework has not yet declared.
#'
#' @param q Numeric scalar or vector in `[0, 1]`. Per-rater accuracy
#'   (`Se = Sp = q`).
#' @param k Integer, odd, `>= 3`. Number of raters voting.
#' @return Numeric vector of EMR_panel values, same length as `q`.
#' @export
#'
#' @examples
#' emr_panel(q = 0.85, k = 5)
#' emr_panel(q = c(0.80, 0.85, 0.90, 0.95), k = 5)
emr_panel <- function(q, k) {
  if (!is.numeric(q) || any(!is.finite(q))) {
    stop("`q` must be a finite numeric vector.", call. = FALSE)
  }
  if (any(q < 0 | q > 1)) {
    stop("`q` must lie in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 3 || k %% 2 == 0 || k != as.integer(k)) {
    stop("`k` must be a single odd integer >= 3.", call. = FALSE)
  }
  k <- as.integer(k)
  threshold <- (k + 1L) %/% 2L   # ceiling(k/2) for odd k
  # P(panel wrong) = P(fewer than `threshold` raters are correct).
  stats::pbinom(threshold - 1L, size = k, prob = q)
}

# ---- classify() ---------------------------------------------------------

#' Column B of the Reporting Card: use-case tier on EMR_panel
#'
#' `classify()` assigns the use-case tier — one of `"Fails"`,
#' `"Indeterminate"`, `"Meets"`, `"Met with distinction"`, or `"Exceeds"` —
#' given the 95% CI on EMR_panel and a declared tolerance `T`. This is
#' Column B of the §8 Reporting Card in the merged GRASS binary-rater-
#' reliability paper.
#'
#' The label sits on EMR_panel (a rate, a function of `q̂`, `F̂`, `k`, and
#' the voting rule) rather than on any individual agreement coefficient.
#' That is what makes the label system regime-conditioned, metric-agnostic,
#' and operationally consequential — the four specific defects of
#' Landis-Koch (1977) and Koo-Li (2016). See `paper2/review/framework_notes.md`
#' §0.6.5 for the full argument.
#'
#' Tier partition (on the 95% CI `[L, U]` of EMR_panel, with tolerance
#' `T`):
#'
#' | Tier | Arithmetic | Meaning |
#' | :--- | :--- | :--- |
#' | `Exceeds` | `U < T_next` | Would pass a more demanding use case |
#' | `Met with distinction` | `U < T/2` (and not Exceeds) | Passes with a safety margin within the declared use case |
#' | `Meets` | `T/2 <= U < T` | Passes the declared use case |
#' | `Indeterminate` | `L < T <= U` | CI straddles T; data too thin to call |
#' | `Fails` | `L >= T` | CI entirely at or above T; panel fails the declared use case |
#'
#' If `emr_lower` is `NA` the `Fails` tier collapses into `Indeterminate`
#' (without a lower bound the function cannot distinguish "CI straddles"
#' from "CI entirely above").
#'
#' @param emr_upper Numeric scalar. Upper bound of the 95% CI on EMR_panel.
#' @param emr_lower Optional numeric scalar. Lower bound of the 95% CI on
#'   EMR_panel. `NA` or `NULL` collapses `Fails` into `Indeterminate`.
#' @param use_case One of `"research"`, `"screening"`, `"diagnostic"`,
#'   `"clinical"`, or `"custom"`. Sets the default tolerance `T` from
#'   [grass_use_case_ladder()]; `"custom"` requires `tolerance`.
#' @param tolerance Optional numeric scalar. Overrides the use-case default.
#'   Required when `use_case = "custom"`.
#' @param tolerance_next Optional numeric scalar. Overrides `T_next`
#'   (used by the `Exceeds` tier). Defaults to the next harder use case on
#'   the ladder; `NA` at the top of the ladder disables `Exceeds`.
#' @param emr_point Optional numeric scalar. EMR_panel point estimate.
#'   Stored on the output for the citation sentence; not used in tier
#'   assignment.
#' @param regime Optional `grass_asymmetry` object from [check_asymmetry()].
#'   Attached so that a Reporting Card row carries both Column A and
#'   Column B together.
#'
#' @return An S3 object of class `grass_classification` with fields:
#' - `tier`: character tier label
#' - `tolerance`: `T` used (with source)
#' - `tolerance_next`: `T_next` used (or `NA`)
#' - `emr_upper`, `emr_lower`, `emr_point`: inputs, echoed for reporting
#' - `margin`: `T - emr_upper` (positive means the panel sits below T)
#' - `use_case`: the use-case label
#' - `regime`: attached Column A (or `NULL`)
#'
#' @seealso [check_asymmetry()] for Column A; [emr_panel()] to compute the
#'   EMR_panel point estimate from `q̂` and `k` under the Tier 1 diagonal.
#' @export
#'
#' @examples
#' # A screening panel, symmetric raters, tight CI
#' asym <- check_asymmetry(se = c(0.86, 0.88, 0.84, 0.87, 0.85),
#'                         sp = c(0.85, 0.87, 0.86, 0.86, 0.87))
#' classify(emr_upper = 0.07, emr_lower = 0.05, use_case = "screening",
#'          emr_point = 0.06, regime = asym)
#'
#' # Same panel held to the diagnostic tolerance — indeterminate
#' classify(emr_upper = 0.07, emr_lower = 0.03, use_case = "diagnostic",
#'          emr_point = 0.05)
#'
#' # Custom tolerance
#' classify(emr_upper = 0.04, use_case = "custom", tolerance = 0.08)
classify <- function(emr_upper,
                     emr_lower = NA_real_,
                     use_case = c("research", "screening", "diagnostic",
                                  "clinical", "custom"),
                     tolerance = NULL,
                     tolerance_next = NULL,
                     emr_point = NA_real_,
                     regime = NULL) {
  use_case <- match.arg(use_case)

  if (!is.numeric(emr_upper) || length(emr_upper) != 1L || !is.finite(emr_upper)) {
    stop("`emr_upper` must be a single finite numeric value.", call. = FALSE)
  }
  if (emr_upper < 0 || emr_upper > 1) {
    stop("`emr_upper` must lie in [0, 1].", call. = FALSE)
  }
  if (!is.null(emr_lower) && length(emr_lower) == 1L &&
      !is.na(emr_lower) && is.finite(emr_lower)) {
    if (emr_lower < 0 || emr_lower > emr_upper) {
      stop("`emr_lower` must lie in [0, emr_upper].", call. = FALSE)
    }
    L <- as.numeric(emr_lower)
  } else {
    L <- NA_real_
  }

  if (use_case == "custom") {
    if (is.null(tolerance)) {
      stop("`use_case = \"custom\"` requires an explicit `tolerance`.",
           call. = FALSE)
    }
    T_val <- as.numeric(tolerance)
    T_source <- "custom (user-declared)"
  } else {
    ladder <- .grass_use_case_ladder[[use_case]]
    T_val <- if (!is.null(tolerance)) as.numeric(tolerance) else ladder$T
    T_source <- if (!is.null(tolerance)) {
      paste0("custom override of ", use_case, " default (", ladder$T, ")")
    } else {
      ladder$source
    }
  }
  if (!is.finite(T_val) || T_val <= 0 || T_val > 1) {
    stop("`tolerance` must be in (0, 1].", call. = FALSE)
  }

  T_next <- if (!is.null(tolerance_next)) {
    as.numeric(tolerance_next)
  } else if (use_case == "custom") {
    NA_real_
  } else {
    .lookup_tolerance_next(use_case)
  }
  if (!is.na(T_next) && (T_next <= 0 || T_next >= T_val)) {
    stop("`tolerance_next` must lie in (0, tolerance).", call. = FALSE)
  }

  U <- as.numeric(emr_upper)
  half_T <- T_val / 2

  # Tier partition. Check from strictest (best) tier down; `Exceeds` wins
  # over `Met with distinction` when both apply (T_next can be above or
  # below T/2 depending on the ladder step).
  tier <- if (!is.na(T_next) && U < T_next) {
    "Exceeds"
  } else if (U < half_T) {
    "Met with distinction"
  } else if (U < T_val) {
    "Meets"
  } else if (is.na(L) || L < T_val) {
    "Indeterminate"
  } else {
    "Fails"
  }

  margin <- T_val - U

  if (!is.null(regime) && !inherits(regime, "grass_asymmetry")) {
    stop("`regime` must be a `grass_asymmetry` object (from `check_asymmetry()`).",
         call. = FALSE)
  }

  out <- list(
    tier = tier,
    use_case = use_case,
    tolerance = T_val,
    tolerance_source = T_source,
    tolerance_next = T_next,
    emr_upper = U,
    emr_lower = L,
    emr_point = if (is.numeric(emr_point) && length(emr_point) == 1L &&
                    is.finite(emr_point)) as.numeric(emr_point) else NA_real_,
    margin = margin,
    regime = regime
  )
  class(out) <- c("grass_classification", "list")
  out
}

#' @export
print.grass_classification <- function(x, digits = 3, ...) {
  cat("grass classification (Column B of Reporting Card)\n",
      "  use case             : ", x$use_case, "\n",
      "  tolerance T          : ", format(x$tolerance, nsmall = 2),
      "   (", x$tolerance_source, ")\n", sep = "")
  if (!is.na(x$tolerance_next)) {
    cat("  T_next (harder use)  : ", format(x$tolerance_next, nsmall = 2), "\n",
        sep = "")
  } else {
    cat("  T_next (harder use)  : n/a (top of ladder or custom)\n", sep = "")
  }
  if (!is.na(x$emr_point)) {
    cat(sprintf("  EMR_panel (point)    : %.*f\n", digits, x$emr_point))
  }
  emr_ci <- sprintf("[%s, %.*f]",
                    if (is.na(x$emr_lower)) "NA" else sprintf("%.*f", digits, x$emr_lower),
                    digits, x$emr_upper)
  cat("  EMR_panel 95% CI     : ", emr_ci, "\n",
      sprintf("  margin T - U         : %+.*f\n", digits, x$margin),
      sep = "")
  cat("  tier                 : ", x$tier, "\n", sep = "")
  if (!is.null(x$regime)) {
    cat("  Column A (model safe): tier ", x$regime$tier,
        "  (", x$regime$regime, ", delta_hat = ",
        format(round(x$regime$delta_hat, digits), nsmall = digits),
        ")\n", sep = "")
  }
  action <- switch(
    x$tier,
    "Fails"                = "Remediate: training, rater reselection, or re-scope use case.",
    "Indeterminate"        = "Collect more data (increase N) or add raters (increase k).",
    "Meets"                = "Proceed.",
    "Met with distinction" = "Proceed; margin is citable evidence of robustness within the declared use case.",
    "Exceeds"              = "Proceed; consider promoting the panel's use case one step up the ladder."
  )
  cat("  action               : ", action, "\n", sep = "")
  invisible(x)
}

#' Coerce a `grass_classification` result to a one-row data.frame
#'
#' @param x A `grass_classification` object.
#' @param row.names,optional,... Standard arguments; ignored except
#'   `row.names`.
#' @return A one-row data.frame suitable for binding into a Reporting Card
#'   table.
#' @export
as.data.frame.grass_classification <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    use_case = x$use_case,
    tolerance = x$tolerance,
    tolerance_next = x$tolerance_next,
    emr_point = x$emr_point,
    emr_lower = x$emr_lower,
    emr_upper = x$emr_upper,
    margin = x$margin,
    tier = x$tier,
    column_a_tier = if (is.null(x$regime)) NA_integer_ else x$regime$tier,
    column_a_regime = if (is.null(x$regime)) NA_character_ else x$regime$regime,
    row.names = row.names,
    stringsAsFactors = FALSE
  )
}
