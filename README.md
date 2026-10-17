# Coherence in hypothesis testing with adaptive significance levels

This repository contains the R code accompanying the paper

**Coherence in hypothesis testing with adaptive significance levels**

by Alia Garrudo Guirado, Luís Gustavo Esteves, and Victor Fossaluza.

## Contents

The `R` directory contains the code used to reproduce the computational
examples and application presented in the paper:

- `coherence_binomial_examples.R`  
  Binomial calculations for Examples 1 and 4.

- `normal_normal_calibration.R`  
  Normal--normal calculations and figures for Example 3.

- `anorexia_application.R`  
  Analysis of the `MASS::anorexia` data, threshold calibration, rejection
  regions, and the interactive Bayes factor surfaces.

## Data

The anorexia application uses the `anorexia` dataset available in the
`MASS` R package.

The response is defined as `Postwt - Prewt` and the analysis compares the CBT and FT treatment groups.

## Software

The analyses were performed in R.

Required packages include:

- MASS
- ggplot2
- plotly
- htmlwidgets

Each script is self-contained and includes the parameter values used in
the corresponding example or application.

## Reproducibility

Run the scripts in the `R` directory to reproduce the numerical results
and figures reported in the paper.
