# `grass` roadmap

## Where the package is today (v0.1.2)

Binary inter-rater agreement for two raters under a spec-dispatch architecture ready to extend to ordinal / multirater / continuous families. Computes Cohen's κ, PABAK, Gwet's AC1 with analytical prevalence-conditioned reference curves at Se = Sp ∈ {0.70, 0.80, 0.85, 0.90}. GRRAS-compliant manuscript methods paragraph via `grass_methods()`. Cohort-parallel reporting via `grass_report_by(..., .parallel = TRUE)`.

- **212 / 212 tests** pass under `load_all + test_dir`
- **`devtools::check(--no-manual)`**: 0 errors, 0 warnings, 1 note (system clock only)
- **Reviewer state**: Dr. M. Chen and J. Okafor both at **5 / 5** after round 5.
- **Next**: v0.1.3 — Fleiss-style logit CI as default (Chen's carry-forward ask; Wilson-logit coverage dips at p < 0.05, n < 150). v0.2.0 — ordinal / multirater / continuous family implementations, empirical Monte-Carlo bands, `grass_power()`.

## Ecosystem architecture (decided)

**One package. One entry point. Metric families dispatched internally.** The package is `grass`. The sole reporting call is `grass_report(data, spec = grass_spec_binary(), ...)`.

Framework labels (GRASS / TURF / MEADOW / FIELD) are **paper and documentation concepts**, not separate installable packages and not separate function prefixes. A user installs `grass`, calls `grass_report()`, and the package handles whichever metric family the data describes. Future metric families slot into the same entry point via different `spec` constructors.

### The framework taxonomy (for the paper)

| Framework label | Metric family | Status |
|---|---|---|
| **GRASS** | Binary, two raters (κ, PABAK, AC1) | implemented (v0.1.x) |
| **TURF** | Ordinal, two raters (weighted κ, ordinal Krippendorff's α) | planned |
| **MEADOW** | Nominal, more than two raters (Fleiss' κ, nominal Krippendorff's α) | planned |
| **FIELD** | Continuous (ICC families, Lin's CCC) | planned |

These labels show up in the paper (section titles, figure captions, methods prose) and in the package's `?grass_roadmap` help page. They never show up as function prefixes the user has to remember.

### The spec-dispatch contract (all families)

```r
grass_report(data, spec, ...)    # → grass_result
```

- `spec` is an S3 object carrying the metric family and its parameters.
- Today: `grass_spec_binary()` is the only constructor that returns a usable spec. The other three (`grass_spec_ordinal()`, `grass_spec_multirater()`, `grass_spec_continuous()`) return placeholder specs that, when passed to `grass_report()`, error with a pointer to `?grass_roadmap`.
- `grass_result$spec$family` carries the family name, which downstream generics (`print`, `plot`, `as.data.frame`, `tidy`, `grass_methods`) dispatch on.

Shared generics, each a single verb that works across every metric family:

- `grass_methods(result)` — manuscript methods paragraph (the ecosystem's signature feature)
- `grass_format_report(result)` — one-line paper-ready summary
- `plot(result)` — landing plot on the family's reference curve
- `tidy(result)` — long-form data.frame for ggplot / gt

## Round 5 scope (next session)

### Architectural refactor

1. Introduce `grass_spec` class with a `$family` slot and family-specific parameters.
2. `grass_spec_binary(reference_level = 0.85, ...)` — the real constructor.
3. Stub constructors for `grass_spec_ordinal()`, `grass_spec_multirater()`, `grass_spec_continuous()` — return valid placeholder specs that error with a `?grass_roadmap` pointer when passed to `grass_report()`.
4. Refactor internals to route through `compute_agreement(data, spec)` / `reference_for(spec, context)` / `classify_regime(metrics, spec)`. All binary today; extensible tomorrow.
5. Rename internal `thresholds` sysdata to `reference_binary` (preserve loading path).
6. `grass_result$spec` slot added.

### Reviewer 5/5 items

7. **`reference_level = 0.70 / 0.80 / 0.85 / 0.90`** — four Se/Sp bands for the binary reference curve. Implemented via analytical reference (closed-form Youden-J-optimal under conditional independence given Se, Sp, prevalence). Validate against the existing simulation at Se = Sp = 0.85; ship the four-band analytical table as `reference_binary.rda`. Defer full Monte Carlo regeneration to a future round if empirical confidence bands are needed.
8. **`.parallel = TRUE` on `grass_report_by()`** — backed by `future.apply::future_lapply` with a `progressr::progressor` handler. `future` and `progressr` go to Suggests. Function errors cleanly if either is missing.

### "Blow anyone away" feature

9. **`grass_methods(result)`** — shared generic with a binary implementation today. Returns a GRRAS-compliant Markdown methods paragraph for the study, pre-filled with the numbers. When TURF / MEADOW / FIELD families arrive, they implement their own methods paragraph via the same generic; the user's code doesn't change.

### Documentation

10. **`?grass_roadmap`** — help page laying out the framework taxonomy and the current implementation state of each family.
11. **Vignette: "The `grass` ecosystem"** — one page on the one-package-many-families design, what's implemented today, what's planned, and why the API won't change when new families arrive.

### Deferred to v0.2.0+

- Actual TURF / MEADOW / FIELD implementations (ordinal, multi-rater, continuous).
- Full Monte Carlo regeneration of reference tables at multiple quality bands (current round uses analytical curves).
- Simulation infrastructure for non-binary families.

## Open questions for round 5 (no sign-off needed to start)

- Should `grass_methods()` emit Markdown, LaTeX, or plain-text by default? Probably Markdown (copy-paste to an Rmd is the common path). Include a `format = "latex"` option for journals that want LaTeX source.
- Should the spec-dispatch errors for unimplemented families offer a "beta-signup" link or just point at `?grass_roadmap`? Just the help page, keep it clean.

## Reviewer context (carry into next session)

- **Dr. M. Chen (clinical epidemiologist)**: 4.75 / 5. 5/5 condition = `reference_level` with four Se/Sp bands properly baked. Signed off on: `ci_width = FALSE` default (opt-in), rescaled kappa CI cutpoints at 0.10/0.20 half-width, floating-point fix, `?grass_report` sample-size-caveats table.
- **J. Okafor (applied data scientist)**: 4.75 / 5. 5/5 condition = `.parallel = TRUE` on `grass_report_by`. Signed off on: `response =` alias, `tidy()` long-form schema, 19-column `as.data.frame`, `compact = TRUE`, soft `id_col` hint. Round-5 stretch ask: `grass_power()` sample-size helper.

Both reviewers happy to continue reviewing through the middleman R-execution workflow when new probes are needed.
