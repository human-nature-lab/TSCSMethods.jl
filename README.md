# TSCSMethods.jl

[![CI](https://github.com/human-nature-lab/TSCSMethods.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/human-nature-lab/TSCSMethods.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/human-nature-lab/TSCSMethods.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/human-nature-lab/TSCSMethods.jl/actions/workflows/Documentation.yml)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://human-nature-lab.github.io/TSCSMethods.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://human-nature-lab.github.io/TSCSMethods.jl/dev)

**Matching methods for causal inference with time-series cross-sectional data**

TSCSMethods.jl implements non-parametric generalized difference-in-differences estimation with covariate matching for panel data. The package provides tools for causal inference in staggered treatment designs, where units receive treatment at different times.

## Key Features

- **Staggered treatment designs**: Handle units treated at different times
- **Covariate matching**: Match treated units to similar controls using time-varying covariates  
- **Flexible time windows**: Specify pre-treatment matching periods and post-treatment estimation periods
- **Multiple balancing strategies**: Manual and automatic covariate balancing
- **Bootstrap inference**: Weighted block-bootstrap for uncertainty quantification
- **Extensions**: Calipers, stratification, refinement, spillover effects
- **Event studies**: Focus on treatment effects relative to event timing

## Quick Start

```julia
using TSCSMethods

# Load example data
data = example_data()

# Create model: match on 10 pre-treatment periods, estimate 5 post-treatment effects
model = makemodel(data, :t, :id, :gub, :Y, [:X1, :X2], -10:-1, 1:5)

# Perform matching and balancing
match!(model)
autobalance(model)

# Estimate treatment effects
estimate!(model, dobayesfactor=false)

# View results
model.resultsoverall.att  # Average treatment effects
```

## Installation

TSCSMethods.jl requires Julia 1.9 or later. Install from the Julia REPL:

```julia
using Pkg
Pkg.add("https://github.com/human-nature-lab/TSCSMethods.jl")
```

## Documentation

- [**Tutorial**](https://human-nature-lab.github.io/TSCSMethods.jl/dev/tutorial/): Step-by-step analysis walkthrough
- [**Methodology**](https://human-nature-lab.github.io/TSCSMethods.jl/dev/methodology/): Statistical methods and assumptions  
- [**API Reference**](https://human-nature-lab.github.io/TSCSMethods.jl/dev/api/): Complete function documentation

## Examples

See the [Jupyter notebook vignette](./vignette/vignette.ipynb) for a complete analysis example with simulated data.

## Method Overview

The package implements the matching approach from Imai et al. (2021) for time-series cross-sectional data:

1. **Matching**: For each treated unit, find control units with similar covariate histories
2. **Balancing**: Assess and improve covariate balance between treated and control groups  
3. **Estimation**: Calculate average treatment effects using matched controls
4. **Inference**: Bootstrap resampling for confidence intervals and significance testing

This approach addresses key challenges in panel data analysis: selection bias, time-varying confounding, and temporal correlation.

## System Requirements

- **Julia**: 1.9 or later
- **Memory**: 8GB+ recommended for large datasets
- **OS**: Windows, macOS, or Linux

The package automatically installs all required dependencies.

## Citation

If you use TSCSMethods.jl in your research, please cite:

```bibtex
@article{feltham_mass_2023,
  title={Mass gatherings for political expression had no discernible association with the local course of the COVID-19 pandemic in the USA in 2020 and 2021},
  author={Feltham, Eric and Forastiere, Laura and Alexander, Marcus and Christakis, Nicholas A},
  journal={Nature Human Behaviour},
  year={2023},
  publisher={Nature Publishing Group}
}

@misc{feltham_tscsmethods_2023,
  title={TSCSMethods.jl: Matching methods for causal inference with time-series cross-sectional data},
  author={Feltham, Eric Martin},
  year={2023},
  url={https://github.com/human-nature-lab/TSCSMethods.jl}
}
```

Please also cite the foundational methodology:

```bibtex
@article{imai_matching_2021,
  title={Matching Methods for Causal Inference with Time-Series Cross-Sectional Data},
  author={Imai, Kosuke and Kim, In Song and Wang, Erik H},
  journal={American Journal of Political Science},
  year={2021},
  publisher={Wiley Online Library}
}
```

## Contributing

Contributions are welcome! Please see the [documentation](https://human-nature-lab.github.io/TSCSMethods.jl/dev/) for development guidelines.

## References

- Imai, K., Kim, I. S., & Wang, E. H. (2021). Matching Methods for Causal Inference with Time-Series Cross-Sectional Data. *American Journal of Political Science*.
- Feltham, E., Forastiere, L., Alexander, M., & Christakis, N. A. (2023). Mass gatherings for political expression had no discernible association with the local course of the COVID-19 pandemic in the USA in 2020 and 2021. *Nature Human Behaviour*.
- Kim, I. S., Ruah, A., Wang, E., & Imai, K. (2020). Insongkim/PanelMatch [R, C]. https://github.com/insongkim/PanelMatch (Original work published 2018)