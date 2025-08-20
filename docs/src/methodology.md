# Methodology

This page explains the statistical methodology implemented in TSCSMethods.jl.

## Overview

TSCSMethods implements the matching approach for time-series cross-sectional (TSCS) data developed by Imai et al. (2021). This method addresses key challenges in causal inference with panel data:

1. **Selection bias**: Units self-select into treatment
2. **Time-varying confounding**: Confounders change over time
3. **Temporal correlation**: Outcomes are correlated within units over time

## The Matching Framework

### Problem Setup

Consider panel data with:
- **Units**: $i = 1, ..., N$ (e.g., counties, countries)
- **Time periods**: $t = 1, ..., T$
- **Treatment**: $D_{it} \in \{0, 1\}$ 
- **Outcome**: $Y_{it}$
- **Covariates**: $X_{it}$

### Staggered Treatment Design

TSCSMethods handles **staggered adoption** where units receive treatment at different times:
- Unit $i$ receives treatment at time $\tau_i$
- $D_{it} = 1$ if $t = \tau_i$, $0$ otherwise
- Focus on **event studies**: effects relative to treatment timing

### Matching Algorithm

For each treated unit $i$ at time $\tau_i$:

1. **Define matching window**: Pre-treatment periods $L = \{l_1, ..., l_L\}$ (negative values)
2. **Find similar controls**: Units $j$ with similar $X_{j,\tau_i+l}$ for $l \in L$
3. **Calculate distances**: $d(i,j) = \sum_{l \in L} w_l ||X_{i,\tau_i+l} - X_{j,\tau_i+l}||$
4. **Select matches**: Closest $K$ control units for each treated unit

### Balancing

Assess matching quality by comparing covariate distributions:

$$\\text{Balance}_{l,k} = \\frac{1}{n_1} \\sum_{i \\in \\text{Treated}} X_{i,\\tau_i+l} - \\frac{1}{n_0} \\sum_{j \\in \\text{Controls}} X_{j,\\tau_i+l}$$

Good matches have small balance statistics.

### Treatment Effect Estimation

Estimate Average Treatment Effect on Treated (ATT) for post-treatment periods $f \\in F$:

$$\\hat{\\tau}_f = \\frac{1}{|\\text{Treated}|} \\sum_{i \\in \\text{Treated}} \\left[ Y_{i,\\tau_i+f} - \\frac{1}{|M_i|} \\sum_{j \\in M_i} Y_{j,\\tau_i+f} \\right]$$

where $M_i$ is the set of matched controls for unit $i$.

### Bootstrap Inference

Uncertainty quantification via bootstrap:

1. **Resample** treated units with replacement
2. **Re-estimate** ATT for each bootstrap sample
3. **Calculate** confidence intervals from bootstrap distribution

## Implementation Details

### Time Period Specification

- **L periods**: Pre-treatment periods for matching (negative values)
  - Example: `L = -10:-1` uses 10 periods before treatment
- **F periods**: Post-treatment periods for estimation (positive values)  
  - Example: `F = 1:5` estimates effects 1-5 periods after treatment
- **Reference period**: Usually `-1` (period just before treatment)

### Covariate Handling

- **Time-invariant**: Covariates constant over time (e.g., geography)
- **Time-varying**: Covariates that change (e.g., population, economic indicators)
- **Standardization**: Covariates standardized by treated unit standard deviation

### Distance Metrics

Default uses weighted Euclidean distance with standardization:
$$d(i,j) = \\sum_{l \\in L} \\sum_{k} \\frac{(X_{ikl} - X_{jkl})^2}{\\sigma^2_{k,l}}$$

where $\\sigma^2_{k,l}$ is the variance of covariate $k$ in period $l$ among treated units.

## Assumptions

The method relies on several key assumptions:

1. **Unconfoundedness**: $Y_{it}(0), Y_{it}(1) \\perp D_{it} | X_{it}, \\text{past}$
2. **Common support**: Sufficient overlap in covariate distributions
3. **No anticipation**: Units don't change behavior before treatment
4. **SUTVA**: No spillover effects between units

## Extensions

TSCSMethods supports several extensions:

- **Calipers**: Restrict matches to units within distance threshold
- **Stratification**: Separate analysis by subgroups
- **Multiple outcomes**: Analyze several dependent variables
- **Refinement**: Iterative improvement of matches

## References

- Imai, K., Kim, I. S., & Wang, E. H. (2021). Matching Methods for Causal Inference with Time-Series Cross-Sectional Data. *American Journal of Political Science*.
- Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity score in observational studies for causal effects. *Biometrika*, 70(1), 41-55.