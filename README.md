# TSCSMethods.jl


[![CI](https://github.com/emfeltham/TSCSMethods.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/emfeltham/TSCSMethods.jl/actions/workflows/CI.yml)
[![Format Check](https://github.com/emfeltham/TSCSMethods.jl/actions/workflows/Format.yml/badge.svg)](https://github.com/emfeltham/TSCSMethods.jl/actions/workflows/Format.yml)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://human-nature-lab.github.io/TSCSMethods.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://human-nature-lab.github.io/TSCSMethods.jl/dev)
[![Build Status](https://travis-ci.com/emfeltham/TSCSMethods.jl.svg?branch=master)](https://travis-ci.com/emfeltham/TSCSMethods.jl)
[![Coverage](https://codecov.io/gh/emfeltham/TSCSMethods.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/emfeltham/TSCSMethods.jl)

## Contents

- [TSCSMethods](#tscsmethods)
  - [Contents](#contents)
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
  - [Hardware Requirements](#hardware-requirements)
  - [Software Requirements](#software-requirements)
    - [OS Requirements](#os-requirements)
- [Installation Guide](#installation-guide)
    - [Package dependencies](#package-dependencies)
    - [Package Installation](#package-installation)
- [Vignette](#vignette)
- [Citation](#citation)

# Overview

Non-parametric generalized differences-in-differences estimation, with covariate matching.

# Repo Contents

- [src](./src): `julia` package code.
- [docs](./docs): package documentation, and usage of the `TSCSMethods` package.

# System Requirements

## Hardware Requirements

`TSCSMethods` works on a standard computer, with sufficient RAM and processing power to support the size of the dataset analyzed by the user. This will be a computer with at least 16 GB, and 4 cores.

The package was tested on a computer with 64 GB of RAM, 16 cores @ 3.4Ghz.

## Software Requirements

### OS Requirements

This package was tested on on MAC OSX 17.0. All of the underlying dependencies are compatible with Windows, Mac, and Linux systems.

This package has been tested on Julia 1.7.1.

# Installation Guide

Julia may be installed on Mac OSX using homebrew <https://brew.sh> by executing:

```shell
brew install julia
```

Otherwise, consult the Julia Language website for installation on your system <https://julialang.org/downloads/>.

### Package dependencies

Users should install the following packages prior to installing `TSCSMethods`, from within a `julia` session:

```{julia}
pkgs = ["Random", "DataFrames", "Dates", "CSV", "JLD2"]

import Pkg
for pkg in pkgs; Pkg.add(pkg) end
```

which will install in less than 5 minutes with the recommended specs.

The `TSCSMethods` package functions with all packages in their latest versions as they appear on `CRAN` on March 08, 2022. The versions of all Julia package dependencies (for TSCSMethods) may be found in the "Manifest.toml" file, and are installed with the package automatically.

### Package Installation

From within a `julia` session, type:

```{julia}
import Pkg; Pkg.add("https://github.com/human-nature-lab/TSCSMethods.jl")
```

The package should take approximately 1 minute to install. 

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

