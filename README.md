# grass

**G**uide for **R**ater **A**greement under **S**tructural **S**kew.

An R package for prevalence-aware interpretation of binary inter-rater agreement. Computes Cohen's kappa, PABAK, and Gwet's AC1, and interprets them against simulation-derived thresholds that depend on the observed prevalence of the positive class — rather than the single fixed threshold of Landis & Koch (1977).

## Why

Landis-Koch thresholds were calibrated at balanced prevalence. At skewed prevalence (common in clinical rating tasks, content moderation, rare-event labeling), a fixed κ cutoff misclassifies rater quality in the majority of scenarios. The GRASS framework conditions interpretation on prevalence.

## What this package is building toward

`grass` is the foundation of a single authoritative framework for rater reliability. The binary family (Cohen's κ, PABAK, AC1) is implemented today at v0.1.2. Planned families — ordinal (weighted κ, AC2, α), multi-rater nominal (Fleiss' κ, Krippendorff's α), and continuous (ICC forms, Lin's CCC, Bland-Altman, G-theory) — will dispatch through the same `grass_report()` entry point via different `grass_spec_*()` constructors. **Prevalence-aware reference curves extend to every family**; the framework's unifying commitment is that interpretation of any agreement coefficient must be conditioned on the study's actual parameters, not applied from a fixed 1977 cutoff table. See `?grass_roadmap` and the *grass ecosystem* vignette for the full taxonomy.

## Install (local dev)

```r
# from the package directory:
devtools::load_all()

# or install the built tarball:
# R CMD INSTALL grass_0.1.0.tar.gz
```

## A 30-second example

```r
library(grass)

# 2x2 counts: cell[1,1] = both rate negative, cell[2,2] = both rate positive
tab <- matrix(c(88, 10, 14, 88), nrow = 2,
              dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))

result <- grass_report(tab, format = "matrix")
result
# plot(result)
```

## What `grass_report()` gives you

A contextual profile — metrics in the context of the skew structure that
shaped them:

- Cohen's kappa, PABAK, Gwet's AC1 (plus confidence intervals)
- Prevalence index and bias index (PI, BI) — the two forms of skew
- A regime label (`balanced`, `prevalence-dominated`, `bias-dominated`,
  `mixed`)
- Optional: signed distance from the GRASS reference curve at the
  observed prevalence

Interpretation belongs to the person reading the numbers; the package
supplies the numbers and the context.

`grass_format_report(result)` condenses all of the above into a single
paper-ready line:

```
κ = 0.76 [0.70, 0.82], PABAK = 0.76, AC1 = 0.76, PI = 0.00, BI = 0.02, N = 200, prevalence = 0.50, balanced regime
```

## Input formats

```r
grass_compute(tab, format = "matrix")              # 2x2 count matrix
grass_compute(df, format = "wide")                  # data.frame: one row/subject, two rater columns
grass_compute(df, format = "long")                  # data.frame: subject, rater, rating columns
grass_compute(list(r1, r2), format = "paired")     # two equal-length rating vectors
```

## Status

v0.1.2 (development). Binary family implemented; spec-dispatch architecture in place for future families. See `NEWS.md` for the changelog and `ROADMAP.md` for the development roadmap.
