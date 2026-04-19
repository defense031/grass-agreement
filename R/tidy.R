# broom::tidy methods — registered conditionally in zzz.R so `broom` is not
# a hard dependency.

tidy.grass_metrics <- function(x, ...) {
  v <- x$values
  data.frame(
    term     = names(v),
    estimate = unname(as.numeric(v)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

tidy.grass_reference <- function(x, ...) {
  x$reference
}

tidy.grass_result <- function(x, ...) {
  # 9 rows: 3 metrics x 3 quantities. `value` is the generic column so the
  # `reference` and `distance` rows don't get misread as point estimates
  # with no CI. `conf.low` / `conf.high` populate only when the row is an
  # `estimate` for a metric that has one (kappa).
  v <- x$metrics$values
  d <- x$distance
  ref_get <- function(m, field = "reference") {
    if (is.null(d)) return(NA_real_)
    out <- d[[field]][d$metric == m]
    if (length(out) == 0) NA_real_ else unname(out)
  }
  est <- c(kappa = unname(v["kappa"]),
           PABAK = unname(v["PABAK"]),
           AC1   = unname(v["AC1"]))
  ref <- c(kappa = ref_get("kappa"),
           PABAK = ref_get("PABAK"),
           AC1   = ref_get("AC1"))
  dist <- c(kappa = ref_get("kappa", "distance"),
            PABAK = ref_get("PABAK", "distance"),
            AC1   = ref_get("AC1", "distance"))

  lo <- c(kappa = unname(v["kappa_wilson_lower"]),
          PABAK = NA_real_,
          AC1   = NA_real_)
  hi <- c(kappa = unname(v["kappa_wilson_upper"]),
          PABAK = NA_real_,
          AC1   = NA_real_)

  metrics <- c("kappa", "PABAK", "AC1")
  data.frame(
    metric    = rep(metrics, times = 3),
    quantity  = rep(c("estimate", "reference", "distance"), each = 3),
    value     = c(est,     ref,      dist),
    conf.low  = c(lo,       rep(NA_real_, 3), rep(NA_real_, 3)),
    conf.high = c(hi,       rep(NA_real_, 3), rep(NA_real_, 3)),
    n         = x$metrics$n,
    prevalence = x$prevalence,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
