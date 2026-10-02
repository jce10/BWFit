Permutation 2: 2fixedC_freeGD_fixedGaussians
============================================

Systematic question
-------------------
What happens to the fitted widths and interference phase when the trapped and
superradiant centroids are prevented from moving?

Constraints
-----------
FIXED:
  - all Gaussian centroids and widths
  - high-angle contaminant centroids and widths
  - trapped centroid E_T at its ground-zero best-fit value
  - superradiant centroid E_SR at its ground-zero best-fit value
  - BG_E0 = 8 MeV

FREE:
  - all component FAC/yield parameters
  - trapped width Gamma_T
  - superradiant width Gamma_SR
  - interference phase delta
  - BG_A0 and BG_A1

Initialization
--------------
Every free parameter is seeded from the corresponding ground-zero best-fit
solution. No calibration, fit interval, binning, background form, or model
component is changed from Permutation 1.

Interpretive goal
-----------------
Permutation 1 showed strong correlations between delta and the Breit-Wigner
centroids, especially rho(delta,E_T)>0 and, at backward angles,
rho(delta,E_SR)<0. This permutation removes both centroid directions. The key
question is whether the angle dependence of delta survives and which remaining
parameters absorb the lost centroid freedom.
