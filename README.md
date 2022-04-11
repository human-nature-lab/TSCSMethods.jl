# TSCSMethods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://emfeltham.github.io/tscsmethods.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://emfeltham.github.io/tscsmethods.jl/dev)
[![Build Status](https://travis-ci.com/emfeltham/tscsmethods.jl.svg?branch=master)](https://travis-ci.com/emfeltham/tscsmethods.jl)
[![Coverage](https://codecov.io/gh/emfeltham/tscsmethods.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/emfeltham/tscsmethods.jl)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Vignette](#demo)
- [Citation](#citation)

# Overview

Non-parametric generalized differences-in-differences estimation, with covariate matching.

# Repo Contents

- [src](./src): `julia` package code.
- [docs](./docs): package documentation, and usage of the `TSCSMethods` package.
- [test](./tests): `julia` unit tests.

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

# Vignette

### execution

After installation of the software dependencies, a given model should take anywhere from 10-60 minutes to execute, depending on workstation specifications.

The models for a given scenario may be executed together, by executing the "run" file, e.g., "run_blm.jl", or separately, by executing the file for a specific model (e.g., "base_model.jl"). Each scenario is housed in a specific subdirectory.

The output for each model is saved as a Julia Data Format object
(https://juliaio.github.io/JLD2.jl/dev/), and the file structure itself
depends on the above packages.

Executing a "run" file creates output files for each model in the "out", subdirectories for each scenario.

You may inspect model output according to the following:

```{julia}
using TSCSMethods, JLD2

output = load_object("ga full_death_rte_.jld2")
```

the output object contains the following fields:
* model
* refinedmodel
* calmodel
* refcalmodel
* matchinfo
* obsinfo

For example, access the results of the refined caliper model as:

```{julia}
output.refcalmodel.results
```

Which yields the ATT estimates, confidence intervals, and the
number of treated units from a given model.

The latter two fields contain information about the matched units and
the treated observation units.

For a simple program example, to estimate the ATTs for a specific event, in
a simple context, see "ga-election/base_model.jl" which runs through
estimation of the overall ATTs for the Georgia special election, each scenario
has the same overall structure:

1. Load the packages, data, and set the parameters for estimation:

```{julia}
# preamble.jl

using Random, TSCSMethods, COVIDPoliticalEvents, Dates
import JLD2

Random.seed!(2019)

savepath = "out/";
datapath = "data/";
scenario = "example "

refinementnum = 5; iters = 500;

treatment = :event; outcome = :outcome;

# clean data with integer-valued binary treatment variable
dat = JLD2.load_object(datafile);

matchingcovariates = [:covar1, :covar2]
timevary = Dict(:covar1 => true, :covar2 => false)

makemodel(
  dat, t, id, treatment, outcome,
  matchingcovariates, timevary,
  10:20, -10:1:-1;
  title = "model",
  estimator = "ATT"
)

dat = dataprep(dat, model);
```

2. Perform matching:

```{julia}
match!(model, dat);
```

3. Perform balancing:

```{julia}
balance!(model, dat);
```

4. Estimate the ATTs for the non-refined, non-caliper model:

```{julia}
estimate!(model, dat);
```

5. Refine the model to the `refinementnum` (5) best matches, and estimate:

```{julia}
refinedmodel = refine(
  model, dat;
  refinementnum = refinementnum, dobalance = true, doestimate = true
);
```

6. Specify an initial caliper:

```{julia}
ibs = Dict(:covariate1 => 0.25)
```

7. Successively apply calipers to get a refined caliper model with sufficient balance:

```
calmodel, refcalmodel, overall = autobalance(
  model, dat;
  calmin = 0.08, step = 0.05,
  initial_bals = ibs,
  dooverall = true
);

8. Save a record of the model:

```{julia}
recordset = makerecords(
  dat, savepath, [model, refinedmodel, calmodel, refcalmodel]
)

# save the average of the ATT values, over the F range (10:20)
JLD2.save_object(savepath * "overall_estimate.jld2", overall)
```

# Citation

For usage of the package and associated manuscript, please cite as:

TSCSMethods
