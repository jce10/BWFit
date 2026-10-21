# Common-reference BW systematics, part 2

45 input cards: 10, 15, 20, 45, and 50 degrees, with nine configurations.
All paths are relative to the BWFit repository root. Copy/merge this bundle's
config/, outputs/, and results/ directories into your repository.

## Run one angle or one permutation

From the BWFit root, run P1 at 10 degrees:

```bash
bash config/thesis/dp_systematics/sys_study2_commonsetpoint/1freeCGD_fixedGaussians_pt2/run_1freeCGD_fixedGaussians_pt2.sh ./bin/bwgaussfit

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/2fixedC_freeGD_fixedGaussians_pt2/run_2fixedC_freeGD_fixedGaussians_pt2.sh

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/3freeC_fixedG_freeD_fixedGaussians_pt2/run_3freeC_fixedG_freeD_fixedGaussians_pt2.sh

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/4fixedCG_freeD_fixedGaussians_pt2/run_4fixedCG_freeD_fixedGaussians_pt2.sh

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/5freeET_fixedGT_ES_GS_freeD_fixedGaussians_pt2/run_5freeET_fixedGT_ES_GS_freeD_fixedGaussians_pt2.sh

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/6fixedET_freeGT_fixedES_GS_freeD_fixedGaussians_pt2/run_6fixedET_freeGT_fixedES_GS_freeD_fixedGaussians_pt2.sh

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/7fixedET_GT_freeES_fixedGS_freeD_fixedGaussians_pt2/run_7fixedET_GT_freeES_fixedGS_freeD_fixedGaussians_pt2.sh

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/8fixedET_GT_ES_freeGS_freeD_fixedGaussians_pt2/run_8fixedET_GT_ES_freeGS_freeD_fixedGaussians_pt2.sh

bash config/thesis/dp_systematics/sys_study2_commonsetpoint/9freeCG_fixedD_fixedGaussians_pt2/run_9freeCG_fixedD_fixedGaussians_pt2.sh
```

Omit the final `10` to run all five angles in that permutation. Each of the
nine permutation directories has its own analogous runner. The optional first
argument selects a different fitter executable. Runners answer the fitter's
JSON prompt automatically and create the ROOT/JSON output directories.
Re-running the same card uses the same output paths and overwrites its results.

For a manual run, call `./bin/bwgaussfit` with the card path, answer `y`, and
paste the result_json path from card_output_index_pt2.csv at the filename prompt.
The card sets the ROOT output; JSON output is supplied by the runner, because
the current reader has no JSON filename field.

## Shared reference

State3: FAC=100, M=7.650 MeV, G=0.060 MeV.
State4: FAC=2000, M=8.200 MeV, G=1.100 MeV.
FAC remains free. All 45 fits start with DELTA=0 in pi/2 units.
P9 fixes DELTA=0 with FIT_PHASE=0 and lets both BW centroids/widths float.
Fixing the phase preserves interference (TOG_INTERFERENCE=1).
The zero fixed phase is the package's explicit reference choice; it was not
specified numerically by the user and is not taken from any fitted phase.

## Fixed/free matrix

| P | E_T | G_T | E_SR | G_SR | delta |
|---|---|---|---|---|---|
| 1 | free | free | free | free | free |
| 2 | fixed | free | fixed | free | free |
| 3 | free | fixed | free | fixed | free |
| 4 | fixed | fixed | fixed | fixed | free |
| 5 | free | fixed | fixed | fixed | free |
| 6 | fixed | free | fixed | fixed | free |
| 7 | fixed | fixed | free | fixed | free |
| 8 | fixed | fixed | fixed | free | free |
| 9 | free | free | free | free | fixed |

Gaussian M/G are fixed to each angle's ground-zero fitted values; Gaussian
FAC starts at its fitted value and floats. Background A0/A1 starts at the
ground-zero fitted values and floats within the original bounds; E0 stays
fixed. Calibrations, windows, binning, and component counts come directly from
the five supplied JSONs. The 20-degree State2 centroid retains the fitted
7.4770000001 MeV value. M_RANGE entries on fixed Gaussians are retained for
bookkeeping and have no effect while FIT_M=0.

Widths that float retain WIDTH_PERC=1 from the supplied inputs, with bounds
computed by BWFit from the common width seeds. The free BW centroids retain
their angle-specific full-fit-window limits.

## Interpretation and validation

P2 is the common-centroid test; P4 fixes all four intrinsic parameters.
Free parameters in other permutations merely share initial values: they are
not tied across angles. P9 tests this specific zero-phase hypothesis; it is
not a simultaneous fit for an unknown common phase.

Card syntax and all numerical settings were checked using the actual devel
FitConfig.cpp reader at commit 11e43e7cbdbe8cb519f2c2481a15ab0ce054d8cf.
The runner scripts were checked for shell syntax and automatic JSON naming.
These cards have not been fitted here; input spectrum ROOT files are not
included. study_manifest_pt2.json and card_output_index_pt2.csv document every
card, constraint, and expected output. results/common_summary_pt2/ is reserved
for later generated summary tables and figures.
