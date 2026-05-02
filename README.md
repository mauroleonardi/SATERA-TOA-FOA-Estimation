# TOA/FOA MATLAB Project Overview

## Purpose
This folder contains MATLAB scripts and generated outputs for ADS-B TOA/FOA estimation experiments.

Main targets:
- TOA (Time of Arrival) estimation
- FOA (Frequency Offset) estimation
- Monte Carlo performance analysis vs SNR, quantization bits, and receiver bandwidth (via decimation)

The project compares multiple estimators and processing chains under controlled impairments.

Algorithms currently covered:
- FOA estimators:
  - FFT peak-based coarse frequency estimate
  - Classical Rife sub-bin interpolation/refinement
  - Modified Rife (`M_Rifenew.m`) with boundary-aware handling near FFT edges
- TOA estimators:
  - Cross-correlation peak detection
  - Interpolated correlation peak for sub-sample timing
  - Rife-refined correlation peak
  - Threshold/rising-edge methods on the ADS-B envelope

Main simulation settings you can tune:
- SNR sweep range and step
- Quantization resolution (for example 4/8/16-bit cases)
- Decimation factor (effective receiver bandwidth)
- Number of Monte Carlo trials per operating point
- Random TOA offset and random FOA generation ranges
- Noise insertion point (variant-dependent, e.g., before/after decimation)
- Outlier filtering/statistics policy (enabled in the newer Monte Carlo variant)

Typical outputs:
- TOA/FOA error statistics (`mean`, `std`, `rms`) across sweeps
- Surface plots and correlation/FFT figures
- `.mat` result files and optional `.xls` summary tables

## Current folder layout (verified)
Everything is now in the same root folder (`src`) plus one results subfolder:

- MATLAB sources:
  - `main_TOA_FOA_mauro.m`
  - `main_TOA_FOA_mauro3.m`
  - `MontecarloTOAFOA.m`
  - `MontecarloTOAFOA2.m`
  - `MontecarloTOAFOA3.m`
  - `M_Rifenew.m`
  - `deci.m`
  - `deci_mauro.m`
- Generated outputs in `src`:
  - `.fig` plots
  - `.mat` summaries (for example `results_4bit.mat`, `results_8bit.mat`, `results_16bit.mat`, `nuovi_results_16bit.mat`)
  - `.xls` tables (`std*`, `rms*`, `mean*`)
- Additional MAT outputs in `src/results`:
  - `nuovi_results_4bit_100prove.mat`
  - `nuovi_results_8bit_100prove.mat`
  - `nuovi_results_16bit_100prove.mat`

No `Xiao_matlab/` or `tesi_thiao_FOA/Matlab Code/` directories are used in the current workspace structure.

## Processing workflow
1. Build ADS-B payload (DF/CA/address/message + CRC + PPM + preamble)
2. Create high-rate waveform
3. Apply random TOA shift and random FOA
4. Downconvert to I/Q and decimate
5. Add noise (position depends on Monte Carlo variant)
6. Quantize envelope and I/Q
7. Estimate FOA with FFT, Rife, and Modified Rife
8. Estimate TOA with correlation, interpolation, Rife-refined peak, and threshold/rising-edge methods
9. Aggregate `rms`, `std`, and `mean` error metrics

## Script roles and dependencies (verified)

### `main_TOA_FOA_mauro.m`
- Calls `MontecarloTOAFOA2` in SNR/decimation sweeps
- Produces 3D surface plots
- Writes XLS tables (`std*`, `rms*`, `mean*`)
- Saves MAT output with `save("results_4bit")`

### `main_TOA_FOA_mauro3.m`
- Calls `MontecarloTOAFOA3` in SNR/decimation sweeps
- Produces the same type of plots
- XLS writes are currently commented out
- Saves MAT output in subfolder with `save("results\\nuovi_results_4bit_100prove")`

### `MontecarloTOAFOA.m`
- Base Monte Carlo engine
- Computes TOA and FOA metrics with a simpler structure

### `MontecarloTOAFOA2.m`
- Main enhanced Monte Carlo engine used by `main_TOA_FOA_mauro.m`
- Includes TOA/FOA randomization refinements and multiple TOA estimator branches
- Uses `M_Rifenew.m` for modified Rife FOA refinement

### `MontecarloTOAFOA3.m`
- Variant engine used by `main_TOA_FOA_mauro3.m`
- Adds noise after decimation in I/Q branches
- Includes outlier filtering with NaN handling before final stats
- Uses `M_Rifenew.m`

### `M_Rifenew.m`
- Modified Rife estimator with edge handling near FFT boundaries
- Returns sub-bin frequency estimate from sampled signal and `dt`

### `deci.m`, `deci_mauro.m`
- Standalone scripts for quick checks and interactive experiments

## What to run
- Full run with XLS export: `main_TOA_FOA_mauro.m`
- Alternative/newer run (MAT in `results/`, XLS disabled): `main_TOA_FOA_mauro3.m`
- Quick manual tests: `deci.m` or `deci_mauro.m`

## Notes
- Simulations can be heavy for large Monte Carlo counts.
- Generated outputs are intentionally stored in the same folder as scripts (plus `results/`).
- `MontecarloTOAFOA3.m` is called by name from `main_TOA_FOA_mauro3.m`; keep file names and function declarations aligned if you plan future refactors.
