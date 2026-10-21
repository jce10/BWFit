# BW systematics, part 3

45 input cards: nine permutations at 10, 15, 20, 45, and 50 degrees.

Extract this archive into the BWFit repository root, merging config/, outputs/,
and results/. Run commands from that repository root. The new input directory is
config/thesis/dp_systematics/sys_study3_commonsetpoint/.

The nine standalone commands are in run_commands_pt3.txt. Each runs all five
angles for one permutation. To select one angle, append ./bin/bwgaussfit and
the angle, for example:

```bash
bash config/thesis/dp_systematics/sys_study3_commonsetpoint/1freeCGD_fixedGaussians_pt3/run_1freeCGD_fixedGaussians_pt3.sh ./bin/bwgaussfit

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/2fixedC_freeGD_fixedGaussians_pt3/run_2fixedC_freeGD_fixedGaussians_pt3.sh

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/3freeC_fixedG_freeD_fixedGaussians_pt3/run_3freeC_fixedG_freeD_fixedGaussians_pt3.sh

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/4fixedCG_freeD_fixedGaussians_pt3/run_4fixedCG_freeD_fixedGaussians_pt3.sh

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/5freeET_fixedGT_ES_GS_freeD_fixedGaussians_pt3/run_5freeET_fixedGT_ES_GS_freeD_fixedGaussians_pt3.sh

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/6fixedET_freeGT_fixedES_GS_freeD_fixedGaussians_pt3/run_6fixedET_freeGT_fixedES_GS_freeD_fixedGaussians_pt3.sh

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/7fixedET_GT_freeES_fixedGS_freeD_fixedGaussians_pt3/run_7fixedET_GT_freeES_fixedGS_freeD_fixedGaussians_pt3.sh

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/8fixedET_GT_ES_freeGS_freeD_fixedGaussians_pt3/run_8fixedET_GT_ES_freeGS_freeD_fixedGaussians_pt3.sh

bash config/thesis/dp_systematics/sys_study3_commonsetpoint/9freeCG_fixedD_fixedGaussians_pt3/run_9freeCG_fixedD_fixedGaussians_pt3.sh

```

## Common starting values

STATE3: FAC=50, M=7.65 MeV, G=0.050 MeV, IS_BW=1.
STATE4: FAC=1500, M=8.2 MeV, G=1.1 MeV, IS_BW=1.
Both FACs are free fit parameters. FIT_M and FIT_G follow this matrix:

| P | E_T | G_T | E_SR | G_SR | phase |
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

Phase starts at zero in every permutation; P9 fixes it at zero while retaining
interference. Gaussian parameters, backgrounds, calibrations, binning, fit
windows, and all other settings are preserved from the pt2 cards.

WIDTH_PERC=1 is preserved. Consequently a floating trapped width now has bounds
0 to 0.100 MeV; previously it was 0 to 0.120 MeV. In permutations fixing the
trapped width, its value is now 0.050 MeV. This pass changes width constraints
as well as initialization and should be described accordingly.

Names, ROOT files, and JSON results all use pt3. Runners find cards beside
themselves even if their containing directories are moved. Run from the BWFit
root so relative input ROOT paths and output paths resolve correctly.
Repeating a pt3 command overwrites that command's pt3 results.

## Output naming

ROOT: outputs/<permutation_name>/<angle>deg_<permutation_name>.root
JSON: results/<permutation_name>/<angle>deg_<permutation_name>.json

The card sets the ROOT filename. Each runner supplies the matching JSON
filename to the fitter's prompt and creates the output directories.
