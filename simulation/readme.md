# GRASS: Guide for Rater Agreement under Structural Skew

## Overview

GRASS is a four-step interpretation protocol for binary inter-rater agreement metrics. It replaces fixed qualitative labels (e.g., Landis-Koch) with prevalence-conditioned thresholds that account for how prevalence affects Cohen's kappa, PABAK, and Gwet's AC1.

The framework addresses a well-documented problem: the same raters can receive different qualitative labels depending on the prevalence of the condition being rated, even when their diagnostic accuracy is unchanged.

## Simulation

The evidence base is a full-factorial Monte Carlo simulation:

- **38 rater profiles** spanning symmetric, within-rater asymmetric, between-rater asymmetric, and opposite-bias configurations
- **23 prevalence levels** (0.01 to 0.99 at 0.05 spacing, plus 0.48 and 0.52)
- **5 sample sizes** (10, 50, 100, 500, 1000)
- **4,370 scenarios**, each with 5,000–20,000 replications
- **24.1 million total replications**

## Repository Structure

```
R/                          Core simulation functions
  00_packages.R             Package loading
  01_data_generating.R      Binary rating data generator
  02_metrics.R              Agreement metric computation (kappa, PABAK, AC1, etc.)
  03_ground_truth.R         Theoretical ground truth and quality bands
  04_grid.R                 Parameter grid construction
  05_runner.R               Parallel simulation orchestration
  06_threshold_derivation.R ROC-optimal threshold derivation

scripts/                    Analysis and visualization scripts
  run_unified_sim.R         Main simulation (config_unified.yaml)
  run_unified_analysis.R    Threshold derivation, classification accuracy, figures
  run_unified_figures.R     Publication-quality figure generation
  run_unified_reviewer_analyses.R  Jackknife stability, CI coverage, sensitivity

config_unified.yaml         Master configuration (profiles, prevalences, sample sizes)
figures/                    Publication-ready figures
output/                     Simulation results and derived thresholds
```

## Reproducing Results

```bash
# Run the full simulation (~30 minutes on 15 cores)
Rscript scripts/run_unified_sim.R

# Derive thresholds and generate figures
Rscript scripts/run_unified_analysis.R

# Jackknife stability, CI coverage, sensitivity analysis
Rscript scripts/run_unified_reviewer_analyses.R

# Regenerate all publication figures
Rscript scripts/run_unified_figures.R
```

## Key Findings

- Landis-Koch labels produce consistent classifications at balanced prevalence (100% of 1,330 scenarios) but only 36.1% consistency at extreme prevalence (412 of 1,140 scenarios)
- With GRASS prevalence-conditioned thresholds, the three metrics agree on classification in 98.4% of scenarios
- Under bias dominance (BI > PI), all three metrics agree in 100% of scenarios — GRASS does not require bias conditioning
- All 72 discordant scenarios occur in the prevalence-dominant regime at extreme prevalence with small sample sizes

## Authors

Austin D. Semmel and Rachel Gidaro
