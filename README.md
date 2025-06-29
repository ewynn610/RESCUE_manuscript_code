# Simulating Paired and Longitudinal Single-cell RNA Sequencing Data with `rescueSim`

This repository contains code for the analyses in our manuscript:

**Simulating Paired and Longitudinal Single-cell RNA Sequencing Data with rescueSim**  

## Overview

The code in this repository:
- Formats empirical scRNA-seq data
- Simulates paired and longitudinal scRNA-seq data using `rescueSim` and `splatPop`
- Calculates benchmarking metrics (e.g., cluster mixing score, KS statistics)
- Performs power analyses using `rescueSim` simulations  

## Repository Structure

- `scripts/00_*` – Helper functions  
- `scripts/01_*` to `03_*` – Format empirical data  
- `scripts/04_*` to `06_*` – Simulate data (with `rescueSim` and `splatPop`) and calculate benchmarking metrics  
- `scripts/09_*` to `11_*` – Power analysis using `rescueSim` simulations  

## Install `rescueSim`
The `rescueSim` package is available at [https://github.com/ewynn610/rescueSim](https://github.com/ewynn610/rescueSim) and can be installed using:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("https://github.com/ewynn610/rescueSim")
