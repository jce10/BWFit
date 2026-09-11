# BWFit covariance/background patch

Based on the `devel` branch of `jce10/BWFit` as retrieved on 2026-09-11.

Modified files:

- `src/FitModel.cpp`
  - Equal lower/upper bounds are now treated explicitly as fixed parameters with `FixParameter()`.
  - If a configured starting value disagrees with equal bounds (for example `BG_E=6` with `BG_E_MIN=BG_E_MAX=8`), BWFit prints a warning and uses the fixed bound value.
  - The same handling is applied to `BG_E`, `BG_A0`, `BG_A1`, and `BG_A2`.

- `include/FitResults.h`
  - Adds `TFitResult` support and declaration of `PrintFitDiagnostics()`.
  - `WriteFitSummaryJSON()` now accepts the full `TFitResult`.

- `src/FitResults.cpp`
  - Prints fit validity, minimizer status, total/free parameter counts, covariance status, EDM, and chi-square p-value.
  - Prints all free-parameter correlations with the interference phase, sorted by absolute correlation strength.
  - Saves parameter labels, fixed/free flags, full covariance matrix, and full correlation matrix to JSON.
  - Background parameter status is now reported as `FIXED` when equal bounds are used.

- `src/main.cpp`
  - Calls `PrintFitDiagnostics()` after the nominal fit.
  - Passes the retained `TFitResult` into the JSON writer.

Recommended config cleanup:

If `E0` is intended to be fixed at 8 MeV, make the config internally consistent:

```
BG_E:      8.0
BG_E_MIN:  8.0
BG_E_MAX:  8.0
```

The patched code will still catch and correct a mismatch, but keeping the config consistent is cleaner.
