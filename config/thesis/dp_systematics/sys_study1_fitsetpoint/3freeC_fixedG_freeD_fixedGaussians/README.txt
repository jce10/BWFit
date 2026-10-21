Permutation 3: free centroids, fixed widths, free phase
==========================================================

Purpose
-------
Complement Permutation 2 by fixing the trapped and superradiant widths at
their angle-specific ground-zero best-fit values while allowing both
Breit-Wigner centroids and the interference phase to float.

FIXED:
- all Gaussian centroids and widths
- high-angle contaminant centroids and widths
- trapped width Gamma_T
- superradiant width Gamma_SR
- BG_E0

FREE:
- all component FAC/yield parameters
- trapped centroid E_T
- superradiant centroid E_SR
- interference phase delta
- BG_A0 and BG_A1

All free parameters are seeded from the corresponding ground-zero best fit.
No calibration, fit interval, binning, background form, or model component is
changed relative to the previous systematic fits.
