# Breit-Wigner Fitter Refactor: Milestone 1

This is the first modular refactor of the original `bw_singles_fitter.cpp` written by K. Hanselman.

## Build

```bash
make clean
make
```

## Run

```bash
./bin/bwsingf config/angle_20_example.txt
```

## Structure

- `include/LineShapes.h`, `src/LineShapes.cpp`: Gaussian, Breit-Wigner, relativistic Breit-Wigner, quadratic background.
- `include/FitConfig.h`, `src/FitConfig.cpp`: reads labeled config files.
- `include/FitModel.h`, `src/FitModel.cpp`: total fit function, including interference.
- `include/FitResults.h`, `src/FitResults.cpp`: reports fit parameter results and corresponding uncertainties.
- `include/SpectrumIO.h`, `src/SpectrumIO.cpp`: ROOT file/tree loading and linear energy calibration.
- `include/Plotting.h`, `src/Plotting.cpp`: draw histogram, total fit, no-interference fit, background, and components.
- `src/main.cpp`: high-level executable flow.

## Important change from original

The executable accepts a config path:

```bash
./bin/bwsingf config/12Cdp_20deg_lin.txt
```

A sample config file and input ROOT file to 


## How the fitter works
1. Read a user-supplied configuration file.
2. Load a ROOT tree containing focal-plane positions.
3. Convert focal-plane position to excitation energy using a linear calibration.
4. Histogram the excitation spectrum.
5. Construct a model consisting of:
    a. Gaussian and/or Breit-Wigner states,
    b. optional trapped-superradiant interference,
    c. linear or quadratic background.
6. Fit the spectrum using ROOT's Minuit2 minimizer.
7. Extract fitted yields and uncertainties.
8. Produce publication-quality plots and residuals.
