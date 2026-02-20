# Mussel Beds — Phase-Plane Demonstration Script

## Overview
This repository contains a self-contained R script (`script_.R`) that demonstrates
phase-plane dynamics for a simple fitness-population density model. The
script focuses on exploring the vector field, equilibria and stochastic
realizations of a two-variable dynamical system rather than on ingesting large
observational datasets.

## License
- **License**: CC0-1.0
- **Copyright**: 2026 Helmholtz-Zentrum hereon GmbH
- **Contributors**: Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>

## Dependencies
The script uses the following R packages:
- `readxl`
- `quantreg`
- `scales`
- `lubridate`
- `deSolve`
- `GA`
- `latex2exp`

Install missing packages with `install.packages()` before running the script.

## What the script does
- Defines a two-variable dynamical system for fitness `f` and normalized
  population `x` with the functions and parameters below.
- Builds a density-based vector-field visualization using a custom
  `vector_field()` helper.
- Plots the phase plane with equilibrium points and a stability map for
  `r/c` vs `mu/k`.
- Runs time simulations using `deSolve::ode()` (Euler method) with a
  simple stochastic `stress.approx()` driver and plots several perturbed
  realizations to show uncertainty.

## Model and parameters (as in `script_.R`)
- Parameters (defaults in the script): `r = 0.5`, `c = 0.5`, `mu = 0.4`,
  `k = 0.75`, `a = 1`, `b = 0`.
- Dynamics function used in the phase-plane visualisation:

  df/dt = r * f * (1 - f) - c * f

  dx/dt = mu * x * (1 - x) - k * (1 - f)^a * x^b * x

- Equilibrium points computed in the script:

  E1 = (0, 0)
  E2 = (1 - c / r, 0)
  E3 = (0, 1 - k / mu)
  E4 = (1 - c / r, 1 - k / mu * c / r)

## Simulation and plotting notes
- The script generates a 2D density image of the vector field and overlays
  equilibrium points and sample trajectories.
- Time simulations use `t = seq(0, 100, by = 0.1)` and an initial state of
  `f = 0.95`, `x = 0.95`. Several realizations are produced by perturbing the
  parameters with Gaussian noise and integrating with the Euler method via
  `deSolve::ode()`.

## Data
This script is demonstrational and does not rely on external input files by
default; it plots directly to the R graphics device. If you want to adapt the
script to use observational mussel or environmental datasets, update the code
to read from `data/` and add preprocessing steps.

## Usage
Open `script_.R` in R (or RStudio) and run the file. The working directory is
set near the top of the script; change `setwd()` if needed.

```r
source('script_.R')
```

## Files
- `script_.R` — main demonstration script producing the phase-plane plots and
  time simulations.

## Version
- Script Date: 2026

