# grass 0.1.2 (development)

Round 5 — spec-dispatch architecture + reviewer 5/5 items.

## Architecture

* **Spec-dispatch**. `grass_report()` now takes `spec = grass_spec_binary()` as its primary argument. Metric families (binary today; ordinal, multirater, continuous planned) are selected by the spec constructor rather than by family-specific function names. One package, one entry point, extensible to future families without API churn. See `?grass_roadmap`.
* **Stub constructors** `grass_spec_ordinal()`, `grass_spec_multirater()`, `grass_spec_continuous()` return valid placeholder specs so users can write release-ready code; passing one to `grass_report()` errors with a pointer to `?grass_roadmap`.
* **Result objects** now carry a `$spec` slot holding the spec used for the call. All downstream verbs (`print`, `plot`, `as.data.frame`, `tidy`, `grass_methods`, `grass_format_report`) dispatch on `result$spec$family`.

## New features

* **`reference_level` argument** on `grass_spec_binary()` accepts `0.70 / 0.80 / 0.85 / 0.90`. Each band is the analytical Youden-J-optimal expected metric value at Se = Sp = reference_level, in closed form under conditional independence of the two raters given the true class. Previously only `"high"` (= 0.85) and `"medium"` (= 0.70) were available.
* **`grass_methods(result)`** — GRRAS-compliant manuscript methods paragraph pre-filled with study numbers. Supports `format = c("markdown", "latex", "plain")`. Templated so that future metric families slot in without touching the dispatcher.
* **`.parallel = TRUE`** on `grass_report_by()` — backed by `future.apply::future_lapply()` with a `progressr::progressor()` handler. Respects the user's active `future::plan()`. Optional dependencies; clean install-command error if absent.
* **`?grass_roadmap`** help page — framework taxonomy, ecosystem architecture, what is implemented today, what is planned.
* **"The grass ecosystem"** vignette — one page on the one-package-many-families design.

## Breaking / estimand changes

* **Reference curve estimand changed.** Pre-0.1.2 shipped a simulation-derived ROC classification threshold (a classifier cutpoint separating rater scenarios of different ground-truth quality bands). 0.1.2 ships the expected metric value at the Se = Sp calibration point — a cleaner analytical estimand with a closed-form interpretation. Numeric values differ at `reference_level = 0.85` vs pre-0.1.2 `reference = "high"` for most prevalences. Full Monte-Carlo regeneration with empirical confidence bands is deferred to a future release.
* **Legacy `reference = "high" / "medium" / "none"` argument on `grass_report()`** is deprecated (soft — emits a one-time message). `"high"` maps to `reference_level = 0.85`, `"medium"` to `0.70`, `"none"` to `NULL`.
* **`grass_reference(quality = ...)`** is deprecated; use `grass_reference(reference_level = ...)`.
* **Internal `thresholds` sysdata renamed** to `reference_binary` (long-form: prevalence × reference_level × metric). `grass_reference_table()` now returns the long form by default and accepts an optional `reference_level` filter. Package-internal only; user impact limited to anyone previously calling `grass:::thresholds`.

---

# grass 0.1.1

Reviewer-driven polish (round 4).

## New in 0.1.1 (round 4)

* `response =` argument — stats-modelling alias for `rater_cols =` (wide) and `rating =` (long). Conflicting canonical + alias errors with a "disagree" message.
* `grass_format_report(ci_width = TRUE)` — appends `"; kappa CI width = 0.26, wide"` with half-width cutpoints **calibrated for kappa**: `tight < 0.10`, `moderate < 0.20`, `wide >= 0.20`. Default `FALSE` (opt-in). Kappa Wilson-logit CIs run 2-2.5x wider than component Se/Sp CIs at the same N, so these cutpoints are higher than the Se/Sp conventions from diagnostic accuracy reporting.
* `broom::tidy.grass_result()` long-form — 9 rows (3 metrics × 3 quantities: `estimate`, `reference`, `distance`). Columns: `metric`, `quantity`, `value`, `conf.low`, `conf.high`, `n`, `prevalence`. `conf.low`/`conf.high` populated only on kappa estimate. `value` (not `estimate`) is the generic column so reference/distance rows do not read as point estimates without CIs. Plots via `ggplot(aes(metric, value, colour = quantity))` without a reshape.
* `grass_report_by(data, group, ...)` — splits `data` by `group` (bare symbol or string), runs `grass_report()` per subset, returns a tidy data.frame with `.cohort` column carrying the group value. Uses `as.data.frame(compact = TRUE)` internally.
* `as.data.frame(r, compact = TRUE)` argument (already in earlier 0.1.1 draft; surface again here).
* Wide-format error now adds a soft `id_col = ` hint when exactly one column is non-binary: *"Column 'xxx' looks like an identifier; you could try `id_col = 'xxx'`."* No hint when ambiguous (all columns binary-looking).
* `reset_grass_warnings()` — clears the once-per-session message cache so loops over many studies with consistently named columns re-emit the fallthrough warning.
* `?grass_report` now has a "Sample size caveats" section documenting the two N thresholds (`N < 10` warning, `N < 30` print note).

## Breaking

* Default `format` is now `"wide"` (was `"matrix"`) across `grass_report()`, `grass_compute()`, `grass_prevalence()`. Users passing a data.frame no longer need to name the format.
* `grass_format_report()` default is now `ascii = TRUE` (was `FALSE`). Safer for Slack, Markdown, and non-UTF-8 locales. Pass `ascii = FALSE` for a Unicode `kappa` in an RStudio console.
* `validate_prevalence()` now errors on exactly `0` or `1` (degenerate cases with no reference defined). Values in `(0, 0.01)` and `(0.99, 1)` still warn and clamp.

## New

* `id_col = NULL` argument on `format = "wide"` — drop a subject-identifier column before the two-rater detection. Example: `grass_report(df, id_col = "subject_id")`.
* `as.data.frame.grass_result()` now returns 19 columns including `kappa_ref`, `kappa_distance`, `PABAK_ref`, `PABAK_distance`, `AC1_ref`, `AC1_distance`, `reference_quality`, and `regime_note`. Pass `compact = TRUE` to drop `regime_note` when binding many rows.
* `print.grass_result()` prints a "sensitive to sampling noise" caveat in the reference-comparison section when `N < 30`. (This is distinct from the `N < 10` warning on `grass_compute()`, which flags the metrics themselves as unreliable at very small N.)

## Fixes

* When no keyword rule matches, `pick_positive()` now escalates from a `message()` to a one-time `warning()` with explicit override guidance. This prevents silent prevalence inversion when rater levels are clinical terms (e.g., `"abnormal"`, `"normal"`) that fall outside the `yes` / `1` / `true` / `positive` keyword list. Flagged by Dr. M. Chen.
* Three-column wide-format error message now names both fix options (`rater_cols = c(...)` and `id_col = "..."`) and lists the actual column names present, without guessing which are raters. Flagged by J. Okafor.

## Still v0.1.0 entry-points

`grass_compute()`, `grass_report()`, `grass_reference()`, `grass_reference_table()`, `grass_prevalence()`, `grass_format_report()`, `grass_plot()`.


# grass 0.1.0

Initial local draft.

## Public API

* `grass_compute()` — raw metric panel from 2x2 or rater data
* `grass_report()` — contextual profile: metrics + PI + BI + regime + reference
* `grass_reference()` — look up reference curves at a given prevalence
* `grass_reference_table()` — full internal reference-curve table (23 prevalence points, High and Medium calibration levels)
* `grass_prevalence()` — estimate prevalence from rater marginals
* `grass_format_report()` — one-line paper-ready summary of a result
* `grass_plot()` — plot a grass object

## Design

* `grass_report()` returns a contextual profile (metrics + PI + BI + regime + signed distance from the reference curve). Regime is one of `balanced`, `prevalence-dominated`, `bias-dominated`, `mixed`; each carries a short structural-implication note describing what the algebra of that regime forces on the metrics.
* Plotting: `plot.grass_result()` landing plot and regime scatter, `plot.grass_reference()` curves, `theme_grass()`. Inline or legend curve labels via `labels = c("auto", "inline", "legend")`.
* Vignettes: `grass-intro`, `skewed-examples`.

## Notes for early readers

* In the paper, the reference values are called *prevalence-conditioned thresholds* (Youden-J-optimal metric values from the simulation). The package uses the name `reference` throughout its public API to signal that the values serve as comparison points reported alongside PI and BI, rather than as cutoffs for a binary verdict.
