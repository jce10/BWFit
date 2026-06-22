#ifndef BWFIT_UPPERLIMIT_H
#define BWFIT_UPPERLIMIT_H

#include "FitConfig.h"
#include <vector>

class TH1F;

namespace bwfit {

struct UpperLimitPoint {
  double forced_yield = 0.0;
  double chi2 = 0.0;
  int ndf = 0;
  double chi2_ndf = 0.0;
  int fit_status = -999;
};

struct UpperLimitScanResult {
  std::vector<UpperLimitPoint> points;
  double best_chi2 = 0.0;
  int best_ndf = 0;
};

// Scan fixed superradiant yields. For each trial yield, the superrad yield
// parameter is fixed while all other free parameters are refit.
UpperLimitScanResult ScanUpperLimit(
    TH1F& hist,
    const FitConfig& cfg,
    const std::vector<double>& best_params,
    const std::vector<double>& trial_yields,
    bool interference_enabled = false);

// Convenience helper: make a simple linear grid of yields.
std::vector<double> MakeYieldGrid(double ymin, double ymax, int npoints);

// Print a compact scan table to stdout.
void PrintUpperLimitScan(const UpperLimitScanResult& scan);

} // namespace bwfit

#endif // BWFIT_UPPERLIMIT_H
