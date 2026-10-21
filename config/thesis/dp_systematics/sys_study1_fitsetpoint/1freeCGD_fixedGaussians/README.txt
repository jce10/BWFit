Permutation 1: 1freeCGD_fixedGaussians — ground-zero reseeded
=================================================================

Systematic definition
---------------------
Gaussian shape parameters are fixed at their ground-zero fitted values.
Trapped centroid/width, superradiant centroid/width, interference phase,
all state yields, and linear-background A0/A1 remain free.

Numerical seeding convention
----------------------------
Every parameter that remains free is initialized at its corresponding
ground-zero best-fit value:
  - all state FAC values
  - trapped M and G
  - superradiant M and G
  - DELTA (in pi/2 units)
  - BG_A0 and BG_A1

No fit ranges, bounds, calibrations, binning, model components, or fixed/free
choices were changed relative to permutation 1. This reseeding is only to
start MINUIT in the known ground-zero basin.

The 20-degree State2 Gaussian remains fixed at 7.4770000001 MeV, which was
essentially the lower centroid bound of the ground-zero fit. That issue is
reserved for a separate nuisance/calibration systematic.
