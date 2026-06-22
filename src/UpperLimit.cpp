#include "UpperLimit.h"
#include "FitModel.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1F.h>

namespace bwfit {

std::vector<double> MakeYieldGrid(double ymin, double ymax, int npoints) {
  if (npoints <= 0) return {};
  if (npoints == 1) return {ymin};
  if (ymax < ymin) std::swap(ymin, ymax);

  std::vector<double> grid;
  grid.reserve(static_cast<std::size_t>(npoints));

  const double step = (ymax - ymin) / static_cast<double>(npoints - 1);
  for (int i = 0; i < npoints; ++i) {
    grid.push_back(ymin + static_cast<double>(i) * step);
  }
  return grid;
}

UpperLimitScanResult ScanUpperLimit(
    TH1F& hist,
    const FitConfig& cfg,
    const std::vector<double>& best_params,
    const std::vector<double>& trial_yields,
    bool interference_enabled) {

  if (best_params.size() != static_cast<std::size_t>(cfg.NumParams())) {
    throw std::runtime_error(
        "ScanUpperLimit: best_params size does not match cfg.NumParams().");
  }

  const int si = cfg.superrad_index;
  if (si < 0 || si >= static_cast<int>(cfg.states.size())) {
    throw std::runtime_error(
        "ScanUpperLimit: cfg.superrad_index is outside cfg.states.");
  }

  const int yield_index = 3 * si + 0;

  FitModel scan_model(cfg, interference_enabled);
  UpperLimitScanResult scan;
  scan.points.reserve(trial_yields.size());

  for (std::size_t i = 0; i < trial_yields.size(); ++i) {
    const double forced_yield = trial_yields[i];

    // Use a unique name for each TF1 so ROOT does not complain about reused names.
    const std::string fit_name = "upper_limit_fit_" + std::to_string(i);
    TF1 test_fit(fit_name.c_str(), scan_model, cfg.fitmin, cfg.fitmax, cfg.NumParams());

    test_fit.SetParameters(best_params.data());
    test_fit.SetNumberFitPoints(cfg.graph_points);
    test_fit.SetNpx(cfg.graph_points);
    ConfigureParameterLimits(test_fit, cfg);

    // This is the actual injection/forced-yield step:
    // hold the superradiant yield fixed and refit all other free parameters.
    test_fit.FixParameter(yield_index, forced_yield);

    TFitResultPtr result = hist.Fit(&test_fit, "RSQ", "", cfg.fitmin, cfg.fitmax);

    UpperLimitPoint point;
    point.forced_yield = forced_yield;
    point.chi2 = test_fit.GetChisquare();
    point.ndf = test_fit.GetNDF();
    point.chi2_ndf = (point.ndf > 0) ? point.chi2 / static_cast<double>(point.ndf) : 0.0;
    point.fit_status = static_cast<int>(result);

    scan.points.push_back(point);
  }

  if (!scan.points.empty()) {
    auto best = std::min_element(
        scan.points.begin(), scan.points.end(),
        [](const UpperLimitPoint& a, const UpperLimitPoint& b) {
          return a.chi2 < b.chi2;
        });
    scan.best_chi2 = best->chi2;
    scan.best_ndf = best->ndf;
  }

  return scan;
}

void PrintUpperLimitScan(const UpperLimitScanResult& scan) {
  std::cout << "\n=== Superradiant Upper-Limit Scan ===\n";
  std::cout << "Forced yield = fixed superrad yield; all other free parameters refit.\n\n";

  std::cout << std::setw(16) << "Y_forced"
            << std::setw(16) << "chi2"
            << std::setw(10) << "NDF"
            << std::setw(16) << "chi2/NDF"
            << std::setw(16) << "Delta chi2"
            << std::setw(10) << "status"
            << "\n";

  for (const auto& p : scan.points) {
    const double delta_chi2 = p.chi2 - scan.best_chi2;
    std::cout << std::setw(16) << std::fixed << std::setprecision(4) << p.forced_yield
              << std::setw(16) << std::fixed << std::setprecision(4) << p.chi2
              << std::setw(10) << p.ndf
              << std::setw(16) << std::fixed << std::setprecision(4) << p.chi2_ndf
              << std::setw(16) << std::fixed << std::setprecision(4) << delta_chi2
              << std::setw(10) << p.fit_status
              << "\n";
  }

  std::cout << "\nCommon one-parameter thresholds: "
            << "Delta chi2 = 1.00 approx 68%, "
            << "2.71 approx 90% one-sided, "
            << "3.84 approx 95% two-sided.\n";
  std::cout << "Use the scan with residual inspection and background/window variations.\n\n";
}

} // namespace bwfit
