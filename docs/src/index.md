```@meta
CurrentModule = TSCSMethods
```

# TSCSMethods.jl

*Matching methods for causal inference with time-series cross-sectional data*

TSCSMethods.jl implements the matching methodology developed in Feltham et al. (2023), which is based on Imai et al. (2021) for causal inference with time-series cross-sectional (TSCS) data.

This package was initially developed for and used in the analyses of Feltham et al. (2023).

## Quick Start

```julia
using TSCSMethods, DataFrames

# Generate example data
dat = example_data(n_units=50, n_days=60)

# Create model
model = makemodel(
    dat, :day, :fips, :gub, :death_rte,
    [:pop_dens], Dict(:pop_dens => false),
    5:10,    # F: post-treatment periods  
    -15:-10  # L: pre-treatment periods
)

# Run complete workflow
match!(model, dat)      # Find matched control units
balance!(model, dat)    # Calculate covariate balances  
estimate!(model, dat)   # Estimate treatment effects

# View results
model.results
```

## Key Features

- **Staggered Treatment Design**: Handles units treated at different times
- **Matching & Balancing**: Find comparable control units and assess covariate balance
- **Bootstrap Inference**: Robust standard errors and confidence intervals
- **Time-Series Structure**: Explicitly accounts for temporal correlation
- **Multiple Outcomes**: Support for analyzing multiple dependent variables

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/emfeltham/TSCSMethods.jl")
```

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
  url={https://github.com/emfeltham/TSCSMethods.jl}
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

## References

- Imai, K., Kim, I. S., & Wang, E. H. (2021). Matching Methods for Causal Inference with Time-Series Cross-Sectional Data. *American Journal of Political Science*.
- Feltham, E., Forastiere, L., Alexander, M., & Christakis, N. A. (2023). Mass gatherings for political expression had no discernible association with the local course of the COVID-19 pandemic in the USA in 2020 and 2021. *Nature Human Behaviour*.
- Kim, I. S., Ruah, A., Wang, E., & Imai, K. (2020). Insongkim/PanelMatch [R, C]. https://github.com/insongkim/PanelMatch (Original work published 2018)

