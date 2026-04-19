#' The grass ecosystem roadmap
#'
#' `grass` is a single package with a single entry point — [grass_report()] —
#' that dispatches to a metric family based on a `grass_spec` object. Each
#' metric family covers one scale type for inter-rater or test-retest
#' reliability, and each is paired with its own prevalence-conditioned
#' reference curve. Today only the binary family is implemented. The rest
#' of the ecosystem is documented here so users can read the API and
#' understand the broader framework without having to track release notes.
#'
#' @section Framework foundation:
#'
#' A core idea animates every family: **fixed interpretation bands
#' (Landis-Koch 1977 and its descendants) are mathematically invalid
#' across prevalences and across scale types.** The reference curve a
#' result is compared against must be conditioned on the actual study
#' parameters (prevalence, design, scale). `grass` is being built as a
#' single authoritative framework that applies this principle uniformly
#' across scale types — binary, nominal, ordinal, continuous,
#' test-retest, internal consistency. The ecosystem complements rather
#' than replaces existing packages (`psych`, `irrCAC`, `irr`,
#' `SimplyAgree`); where those packages compute coefficients, `grass`
#' attaches a prevalence-aware reference and a methods-paragraph
#' scaffold that carries through to the manuscript.
#'
#' @section Framework taxonomy:
#'
#' Four framework labels appear in the paper and in this documentation.
#' They are organising concepts, not function prefixes. The user-facing
#' API never exposes them — every family is accessed through
#' [grass_report()] with the appropriate spec.
#'
#' \tabular{lll}{
#'   **Label** \tab **Metric family** \tab **Status in 0.1.2** \cr
#'   GRASS \tab Binary, two raters (Cohen's kappa, PABAK, Gwet's AC1). Prevalence-conditioned reference curves at Se = Sp in \{0.70, 0.80, 0.85, 0.90\}. \tab **implemented** — see [grass_spec_binary()] \cr
#'   TURF  \tab Ordinal, two raters (weighted Cohen's kappa with linear and quadratic weights, Gwet's AC2, ordinal Krippendorff's alpha). \tab planned — [grass_spec_ordinal()] placeholder only \cr
#'   MEADOW \tab Nominal, more than two raters (Fleiss' kappa, nominal Krippendorff's alpha). \tab planned — [grass_spec_multirater()] placeholder only \cr
#'   FIELD  \tab Continuous (Shrout-Fleiss ICC forms, Lin's concordance correlation, Bland-Altman limits with parametric and non-parametric bounds, standard error of measurement, generalisability-theory variance components). \tab planned — [grass_spec_continuous()] placeholder only \cr
#' }
#'
#' @section Stable API:
#'
#' The same verbs work across every family. When the planned families
#' ship, existing user code does not change — the only thing that
#' changes is the spec constructor called. This is what "one package,
#' one entry point" buys the user.
#'
#' * [grass_report()] with the appropriate spec — compute and attach
#'   the family's reference.
#' * [grass_methods()] — GRRAS-compliant manuscript methods paragraph
#'   pre-filled with study numbers.
#' * [grass_format_report()] — one-line paper-ready summary.
#' * `plot()` — landing plot on the family's reference curve.
#' * [broom::tidy()] — long-form data.frame of estimates, references,
#'   and distances.
#'
#' @section What is not yet implemented:
#'
#' Until the ordinal, multirater, and continuous families land,
#' constructing one of the placeholder specs
#' ([grass_spec_ordinal()], [grass_spec_multirater()],
#' [grass_spec_continuous()]) is legal so users can write code ready
#' for the release, but passing one to [grass_report()] or
#' [grass_methods()] errors with a pointer back here. The research
#' roadmap for each family includes: non-parametric confidence intervals
#' (logit-transformed Wilson for kappa, Harrell-Davis and
#' Sfakianakis-Verginis quantiles for Bland-Altman), Matthews
#' correlation coefficient and specific-agreement fallbacks in the
#' binary family for rare-event diagnostics, and generalisability
#' theory variance components in the continuous family.
#'
#' @section Paper:
#'
#' The foundational paper for the framework and the binary family is in
#' review (Semmel 202X). FIELD / TURF / MEADOW papers will follow,
#' citing the GRASS paper as the methodological precedent. Each paper
#' accompanies a minor release of this package rather than a new
#' package.
#'
#' @name grass_roadmap
#' @aliases grass-roadmap
NULL
