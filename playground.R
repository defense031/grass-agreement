# grass v0.1.2 playground
# Load the package via load_all() so you can edit + iterate without reinstall.
# Run sections below top-to-bottom, or cherry-pick.

devtools::load_all(".")
library(ggplot2)


# ============================================================
# 1. The new spec-dispatch API
# ============================================================
# grass_report() takes a `spec` argument that selects the metric family
# and its parameters. Today only the binary family is implemented;
# ordinal/multirater/continuous constructors are placeholders that
# point at ?grass_roadmap.

grass_spec_binary()                               # default: reference_level = 0.85
grass_spec_binary(reference_level = 0.90)
grass_spec_binary(reference_level = NULL)         # skip reference attachment

# The stub constructors build cleanly so release-ready code compiles…
grass_spec_ordinal()
# …but passing one to grass_report() errors with a roadmap pointer:
tryCatch(
  grass_report(matrix(c(80, 10, 10, 80), 2),
               format = "matrix",
               spec = grass_spec_ordinal()),
  error = function(e) message("expected: ", conditionMessage(e))
)


# ============================================================
# 2. A small example
# ============================================================
# Cohen 1960 Table 1.
tab <- matrix(c(88, 10, 14, 88), nrow = 2,
              dimnames = list(R1 = c("neg", "pos"),
                              R2 = c("neg", "pos")))

r <- grass_report(tab, format = "matrix")
r                                                 # full printed profile

grass_format_report(r)                            # one-line paper summary
grass_format_report(r, ci_width = TRUE)           # + CI width descriptor


# ============================================================
# 3. Reference bands side-by-side
# ============================================================
# The four analytical bands at the same observed study. Watch the
# `reference` and `delta` columns shift as the rater-quality band moves.

bands <- c(0.70, 0.80, 0.85, 0.90)
side_by_side <- do.call(rbind, lapply(bands, function(q) {
  ri <- grass_report(tab, format = "matrix",
                     spec = grass_spec_binary(reference_level = q))
  data.frame(
    reference_level = q,
    kappa_ref = ri$distance$reference[ri$distance$metric == "kappa"],
    kappa_delta = ri$distance$distance[ri$distance$metric == "kappa"],
    PABAK_ref = ri$distance$reference[ri$distance$metric == "PABAK"],
    AC1_ref   = ri$distance$reference[ri$distance$metric == "AC1"]
  )
}))
print(side_by_side)


# ============================================================
# 4. grass_methods(): manuscript paragraph
# ============================================================
# GRRAS-compliant methods paragraph pre-filled with the study's numbers.
# Three formats: markdown, latex, plain.

cat(grass_methods(r, format = "markdown"), "\n\n")
cat(grass_methods(r, format = "latex"),    "\n\n")
cat(grass_methods(r, format = "plain"),    "\n\n")

# Paste the latex variant directly into a manuscript, or:
# knitr::asis_output(grass_methods(r, format = "markdown"))


# ============================================================
# 5. Plot on the reference curve
# ============================================================
plot(r)                                           # landing plot
plot(r, type = "regime")                          # PI^2 vs BI^2

# At a different band:
r_80 <- grass_report(tab, format = "matrix",
                     spec = grass_spec_binary(reference_level = 0.80))
plot(r_80)


# ============================================================
# 6. Cohort-split reporting, with optional parallel dispatch
# ============================================================
set.seed(42)
site_study <- data.frame(
  xray_id = sprintf("CXR-%03d", 1:120),
  rater_A = sample(c("abnormal", "normal"), 120, replace = TRUE, prob = c(0.3, 0.7)),
  rater_B = sample(c("abnormal", "normal"), 120, replace = TRUE, prob = c(0.3, 0.7)),
  site    = rep(c("A", "B", "C", "D"), each = 30)
)

# Sequential:
grass_report_by(site_study, site,
                id_col = "xray_id", positive = "abnormal")

# Parallel — optional deps (future.apply + progressr). Respects your plan():
# future::plan(future::multisession, workers = 2)
# progressr::with_progress(
#   grass_report_by(site_study, site,
#                   id_col = "xray_id", positive = "abnormal",
#                   .parallel = TRUE)
# )
# future::plan(future::sequential)


# ============================================================
# 7. The skewed regime — where Landis-Koch breaks and GRASS shines
# ============================================================
# Rare-event study: 3% prevalence, two competent raters at Se = Sp ~ 0.88.
set.seed(7)
n <- 300; p <- 0.03; q <- 0.88
truth <- rbinom(n, 1, p)
r1 <- ifelse(truth == 1, rbinom(n, 1, q), rbinom(n, 1, 1 - q))
r2 <- ifelse(truth == 1, rbinom(n, 1, q), rbinom(n, 1, 1 - q))

skew <- grass_report(data.frame(r1, r2),
                     spec = grass_spec_binary(reference_level = 0.85))
skew
# Observe: kappa collapses under the kappa paradox, but PABAK and AC1
# hold their signal. The regime label flags 'prevalence-dominated'.

cat(grass_methods(skew), "\n")


# ============================================================
# 8. Tidy for plotting / further analysis
# ============================================================
# Requires broom installed.
if (requireNamespace("broom", quietly = TRUE)) {
  tidy_r <- broom::tidy(r)
  print(tidy_r)
  ggplot(tidy_r[tidy_r$quantity %in% c("estimate", "reference"), ],
         aes(x = metric, y = value, fill = quantity)) +
    geom_col(position = "dodge") +
    theme_grass() +
    labs(title = "Observed vs reference at q = 0.85")
}


# ============================================================
# 9. The roadmap
# ============================================================
# What's planned:
?grass_roadmap                                    # framework taxonomy
vignette("grass-ecosystem")                       # ecosystem design
vignette("grass-intro")                           # the basics
vignette("skewed-examples")                       # where skew matters
