#' Report agreement metrics in their prevalence and bias context
#'
#' Computes the full metric panel, estimates prevalence, and assembles a
#' contextual profile: observed values and CIs; the prevalence index (PI) and
#' bias index (BI) describing the skew structure of the 2x2 table; the
#' regime classification (`"balanced"`, `"prevalence-dominated"`,
#' `"bias-dominated"`, or `"mixed"`) derived from PI and BI; and the signed
#' distance of each metric from the GRASS reference curve at the observed
#' prevalence.
#'
#' @section Sample size caveats:
#' Two distinct small-N signals, kept separate:
#'
#' | Threshold | Channel                 | Meaning                                                             |
#' | --------- | ----------------------- | ------------------------------------------------------------------- |
#' | `N < 10`  | `warning()` at compute  | The metrics themselves are unreliable; CIs near the bounds.         |
#' | `N < 30`  | print-time note         | Metrics are defensible but reference deltas are sensitive to noise. |
#'
#' A study at N = 20 returns defensible kappa / PABAK / AC1 but prints a
#' caveat next to the reference comparison. A study at N = 8 also gets the
#' compute-level warning on top.
#'
#' @inheritParams grass_compute
#' @param spec A `grass_spec` object selecting the metric family and any
#'   family-specific parameters. Defaults to `grass_spec_binary()`
#'   (reference_level = 0.85). See [grass_spec_binary()] for the
#'   `reference_level` choices, and `?grass_roadmap` for the planned
#'   ordinal / multirater / continuous specs.
#' @param prevalence Optional override for the prevalence used for the
#'   reference curve lookup. If `NULL` (the default), prevalence is estimated
#'   from rater marginals via [grass_prevalence()].
#' @param reference Deprecated in 0.1.2 in favour of the spec. Accepts
#'   `"high"` → `reference_level = 0.85`, `"medium"` → `0.70`, `"none"`
#'   → no reference attached. Emits a one-time soft deprecation warning.
#'
#' @return A `grass_result` S3 object with components:
#'   * `metrics` — a `grass_metrics` object
#'   * `spec` — the `grass_spec` used
#'   * `prevalence`, `prevalence_source`
#'   * `regime` — one of `"balanced"`, `"prevalence-dominated"`,
#'     `"bias-dominated"`, `"mixed"`
#'   * `regime_note` — short description of the regime's structural implication
#'   * `reference` — a `grass_reference` object (or `NULL` if the spec
#'     has `reference_level = NULL`)
#'   * `distance` — data.frame of metric, observed value, reference value,
#'     signed distance (observed - reference).
#' @export
#'
#' @examples
#' tab <- matrix(c(88, 10, 14, 88), nrow = 2,
#'               dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
#' grass_report(tab, format = "matrix")
#'
#' # Select a different reference band via the spec:
#' grass_report(tab, format = "matrix",
#'              spec = grass_spec_binary(reference_level = 0.80))
grass_report <- function(data, format = c("wide", "matrix", "long", "paired"),
                         positive = NULL, prevalence = NULL,
                         spec = grass_spec_binary(),
                         reference = NULL, ...) {
  format <- match.arg(format)
  call   <- match.call()

  # Legacy `reference =` argument. If supplied, it overrides the spec's
  # reference_level and emits a one-time soft deprecation. This preserves
  # existing user code while nudging toward the spec-based API.
  if (!is.null(reference)) {
    if (!identical(reference, "high") && !identical(reference, "medium") &&
        !identical(reference, "none")) {
      stop("`reference` must be one of \"high\", \"medium\", or \"none\". ",
           "Got: ", deparse(reference), ".", call. = FALSE)
    }
    msg_once(paste0("deprecate_reference_", reference),
             paste0("grass: `reference = ", shQuote(reference), "` is deprecated. ",
                    "Use `spec = grass_spec_binary(reference_level = ...)` instead. ",
                    "This call's behaviour is unchanged."))
    if (!inherits(spec, "grass_spec_binary")) {
      stop("Legacy `reference =` argument applies only to the binary family. ",
           "Drop `reference =` and configure the spec directly.", call. = FALSE)
    }
    spec <- grass_spec_binary(reference_level = switch(
      reference,
      high   = 0.85,
      medium = 0.70,
      none   = NULL
    ))
  }

  if (!is_grass_spec(spec)) {
    stop("`spec` must be a grass_spec object. See ?grass_spec_binary.",
         call. = FALSE)
  }

  metrics <- compute_agreement(data, spec, format = format, positive = positive, ...)

  if (is.null(prevalence)) {
    p <- estimate_prevalence(metrics$table)
    prev_source <- "rater_marginals"
  } else {
    p <- validate_prevalence(prevalence)
    prev_source <- "user_supplied"
  }

  reg <- classify_regime(metrics, spec)

  ref  <- reference_for(spec, list(prevalence = p))
  dist <- if (is.null(ref)) NULL else distance_to_reference(metrics, ref)

  new_grass_result(
    metrics = metrics,
    spec = spec,
    prevalence = p,
    prevalence_source = prev_source,
    regime = reg$regime,
    regime_note = reg$note,
    reference = ref,
    distance = dist,
    call = call
  )
}

#' Estimate prevalence of the positive class from rater data
#'
#' Averages the marginal positive rates of the two raters.
#'
#' @inheritParams grass_compute
#'
#' @return A single numeric in `[0, 1]`.
#' @export
#'
#' @examples
#' tab <- matrix(c(88, 10, 14, 88), nrow = 2,
#'               dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
#' grass_prevalence(tab, format = "matrix")
grass_prevalence <- function(data, format = c("wide", "matrix", "long", "paired"),
                             positive = NULL, ...) {
  format <- match.arg(format)
  norm <- normalize_input(data, format = format, positive = positive, ...)
  tab <- if (!is.null(norm$table)) norm$table else build_table(norm$r1, norm$r2)
  estimate_prevalence(tab)
}

# Internal: average marginal positive rate from a 2x2 table.
estimate_prevalence <- function(tab) {
  N <- sum(tab)
  if (N == 0) return(NA_real_)
  p1 <- (tab[2, 1] + tab[2, 2]) / N  # R1 says 1
  p2 <- (tab[1, 2] + tab[2, 2]) / N  # R2 says 1
  (p1 + p2) / 2
}

# Internal: classify the (PI, BI) regime and attach the structural
# implication — what the mathematics forces on the metrics, stated as fact.
#   - balanced              : both indices small
#   - prevalence-dominated  : PI meaningfully exceeds BI
#   - bias-dominated        : BI meaningfully exceeds PI
#   - mixed                 : both indices non-trivial with neither dominant
#
# Called by the S3 method classify_regime.grass_spec_binary() in
# dispatch.R. Named `classify_regime_binary` (not `classify_regime`) so
# the generic of the same name can live in dispatch.R without shadowing.
classify_regime_binary <- function(pi, bi, small = 0.10, gap = 0.10) {
  pi <- unname(pi); bi <- unname(bi)
  if (is.na(pi) || is.na(bi)) {
    return(list(regime = NA_character_, note = NA_character_))
  }
  if (max(pi, bi) < small) {
    return(list(
      regime = "balanced",
      note = paste0(
        "Prevalence and bias indices are both small, so kappa, PABAK, and ",
        "AC1 are algebraically close. Substantial disagreement among the ",
        "three metrics in this regime points to model or data issues rather ",
        "than rater behaviour."
      )
    ))
  }
  if (pi - bi > gap) {
    return(list(
      regime = "prevalence-dominated",
      note = paste0(
        "Expected agreement under marginal independence is high, so kappa ",
        "is algebraically suppressed even when raters are accurate on the ",
        "minority class. PABAK and AC1 use different denominators and are ",
        "less sensitive to this effect; a gap between kappa and PABAK/AC1 ",
        "is a structural property of the prevalence."
      )
    ))
  }
  if (bi - pi > gap) {
    return(list(
      regime = "bias-dominated",
      note = paste0(
        "The two raters endorse the positive class at different marginal ",
        "rates. AC1 conditions on this disagreement and will trail kappa ",
        "and PABAK; the appropriate response is rater calibration, not ",
        "more training on individual cases."
      )
    ))
  }
  list(
    regime = "mixed",
    note = paste0(
      "Prevalence and bias effects are both present. No single metric ",
      "cleanly captures agreement here; report kappa, PABAK, and AC1 ",
      "alongside PI and BI so the reader sees the full structure."
    )
  )
}

# Internal: signed distance from observed metric to the reference curve at the
# evaluation prevalence. Positive = above the reference, negative = below.
distance_to_reference <- function(metrics, reference) {
  v <- metrics$values
  ref <- reference$reference
  observed <- c(kappa = unname(v["kappa"]),
                PABAK = unname(v["PABAK"]),
                AC1   = unname(v["AC1"]))
  data.frame(
    metric    = ref$metric,
    observed  = observed[ref$metric],
    reference = ref$reference,
    distance  = observed[ref$metric] - ref$reference,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
