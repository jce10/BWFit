#include "UpperLimit.h"
#include "FitModel.h"
#include "SpectrumIO.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1F.h>

namespace bwfit {

namespace {

// These values are set by ScanUpperLimit() and then used by
// PrintUpperLimitScan(). Keeping them internal avoids changing UpperLimit.h.
double g_last_bin_width_MeV = std::numeric_limits<double>::quiet_NaN();
double g_last_reference_chi2 = std::numeric_limits<double>::quiet_NaN();
int g_last_reference_ndf = 0;
int g_last_reference_status = -999;

struct CrossingResult {
  bool found = false;
  double yield = 0.0;
};

CrossingResult FindUpperCrossing(const UpperLimitScanResult& scan,
                                 double threshold) {
  CrossingResult result;

  if (scan.points.size() < 2) return result;

  // Locate the minimum chi2 point in the scanned grid. We then search only
  // toward larger forced yields to find the upper confidence-limit crossing.
  const auto min_it = std::min_element(
      scan.points.begin(), scan.points.end(),
      [](const UpperLimitPoint& a, const UpperLimitPoint& b) {
        return a.chi2 < b.chi2;
      });

  const std::size_t min_index =
      static_cast<std::size_t>(std::distance(scan.points.begin(), min_it));

  for (std::size_t i = min_index + 1; i < scan.points.size(); ++i) {
    const auto& left = scan.points[i - 1];
    const auto& right = scan.points[i];

    const double left_delta = left.chi2 - scan.best_chi2;
    const double right_delta = right.chi2 - scan.best_chi2;

    if (left_delta <= threshold && right_delta >= threshold) {
      const double denom = right_delta - left_delta;

      // Linear interpolation between adjacent scan points.
      if (std::fabs(denom) > 0.0) {
        const double fraction = (threshold - left_delta) / denom;
        result.yield =
            left.forced_yield +
            fraction * (right.forced_yield - left.forced_yield);
      } else {
        result.yield = right.forced_yield;
      }

      result.found = true;
      return result;
    }
  }

  return result;
}

} // namespace

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

  g_last_bin_width_MeV = CalculateHistogramBinWidth(cfg);

  // ------------------------------------------------------------------------
  // Establish the reference chi2 with the superradiant yield FREE.
  //
  // Delta chi2 must be measured relative to the unconstrained best fit, not
  // merely relative to the lowest point that happens to lie on the scan grid.
  // ------------------------------------------------------------------------
  TF1 reference_fit(
      "upper_limit_reference_fit",
      scan_model,
      cfg.fitmin,
      cfg.fitmax,
      cfg.NumParams());

  reference_fit.SetParameters(best_params.data());
  reference_fit.SetNumberFitPoints(cfg.graph_points);
  reference_fit.SetNpx(cfg.graph_points);
  ConfigureParameterLimits(reference_fit, cfg);

  TFitResultPtr reference_result =
      hist.Fit(&reference_fit, "RSQ", "", cfg.fitmin, cfg.fitmax);

  g_last_reference_chi2 = reference_fit.GetChisquare();
  g_last_reference_ndf = reference_fit.GetNDF();
  g_last_reference_status = static_cast<int>(reference_result);

  scan.best_chi2 = g_last_reference_chi2;
  scan.best_ndf = g_last_reference_ndf;

  // Use the converged unconstrained solution as the common starting point for
  // every constrained fit.
  std::vector<double> reference_params(cfg.NumParams());
  for (int ipar = 0; ipar < cfg.NumParams(); ++ipar) {
    reference_params[ipar] = reference_fit.GetParameter(ipar);
  }

  for (std::size_t i = 0; i < trial_yields.size(); ++i) {
    const double forced_yield = trial_yields[i];

    const std::string fit_name = "upper_limit_fit_" + std::to_string(i);
    TF1 test_fit(
        fit_name.c_str(),
        scan_model,
        cfg.fitmin,
        cfg.fitmax,
        cfg.NumParams());

    test_fit.SetParameters(reference_params.data());
    test_fit.SetNumberFitPoints(cfg.graph_points);
    test_fit.SetNpx(cfg.graph_points);
    ConfigureParameterLimits(test_fit, cfg);

    // Hold the superradiant normalization fixed and profile all other free
    // parameters.
    test_fit.FixParameter(yield_index, forced_yield);

    TFitResultPtr result =
        hist.Fit(&test_fit, "RSQ", "", cfg.fitmin, cfg.fitmax);

    UpperLimitPoint point;
    point.forced_yield = forced_yield;
    point.chi2 = test_fit.GetChisquare();
    point.ndf = test_fit.GetNDF();
    point.chi2_ndf =
        (point.ndf > 0)
            ? point.chi2 / static_cast<double>(point.ndf)
            : 0.0;
    point.fit_status = static_cast<int>(result);

    scan.points.push_back(point);
  }

  return scan;
}

void PrintUpperLimitScan(const UpperLimitScanResult& scan) {
  std::cout << "\n=== Superradiant Upper-Limit Scan ===\n";
  std::cout
      << "Forced normalization = fixed superrad yield parameter; "
      << "all other free parameters refit.\n";

  if (std::isfinite(g_last_bin_width_MeV) &&
      g_last_bin_width_MeV > 0.0) {
    std::cout << "Bin width: dE = "
              << std::fixed << std::setprecision(6)
              << g_last_bin_width_MeV << " MeV/bin = "
              << 1000.0 * g_last_bin_width_MeV << " keV/bin\n";
    std::cout
        << "Integral counts = forced normalization / dE.\n";
  }

  std::cout << "Reference: unconstrained fit with superrad normalization free.\n"
            << "  chi2 = " << std::fixed << std::setprecision(6)
            << g_last_reference_chi2
            << " | NDF = " << g_last_reference_ndf
            << " | status = " << g_last_reference_status << "\n\n";

  std::cout << std::setw(16) << "Y_norm"
            << std::setw(18) << "Integral_counts"
            << std::setw(16) << "chi2"
            << std::setw(10) << "NDF"
            << std::setw(16) << "chi2/NDF"
            << std::setw(16) << "Delta chi2"
            << std::setw(10) << "status"
            << "\n";

  for (const auto& p : scan.points) {
    const double delta_chi2 = p.chi2 - scan.best_chi2;
    const double integral_counts =
        (std::isfinite(g_last_bin_width_MeV) &&
         g_last_bin_width_MeV > 0.0)
            ? p.forced_yield / g_last_bin_width_MeV
            : std::numeric_limits<double>::quiet_NaN();

    std::cout << std::setw(16)
              << std::fixed << std::setprecision(4)
              << p.forced_yield
              << std::setw(18)
              << std::fixed << std::setprecision(2)
              << integral_counts
              << std::setw(16)
              << std::fixed << std::setprecision(4)
              << p.chi2
              << std::setw(10) << p.ndf
              << std::setw(16)
              << std::fixed << std::setprecision(4)
              << p.chi2_ndf
              << std::setw(16)
              << std::fixed << std::setprecision(4)
              << delta_chi2
              << std::setw(10) << p.fit_status
              << "\n";
  }

  std::cout << "\nInterpolated upper-branch crossings:\n";

  struct Threshold {
    double delta_chi2;
    const char* label;
  };

  const Threshold thresholds[] = {
      {1.00, "68% profile interval"},
      {2.71, "95% one-sided upper limit"},
      {3.84, "95% two-sided interval"}
  };

  for (const auto& threshold : thresholds) {
    const CrossingResult crossing =
        FindUpperCrossing(scan, threshold.delta_chi2);

    std::cout << "  Delta chi2 = "
              << std::fixed << std::setprecision(2)
              << threshold.delta_chi2
              << " (" << threshold.label << "): ";

    if (!crossing.found) {
      std::cout << "not reached within the scan range\n";
      continue;
    }

    const double counts =
        (std::isfinite(g_last_bin_width_MeV) &&
         g_last_bin_width_MeV > 0.0)
            ? crossing.yield / g_last_bin_width_MeV
            : std::numeric_limits<double>::quiet_NaN();

    std::cout << "Y_norm = "
              << std::fixed << std::setprecision(4)
              << crossing.yield
              << " counts*MeV"
              << " | integral = "
              << std::fixed << std::setprecision(1)
              << counts << " counts\n";
  }

  // A lower scan point than the unconstrained reference usually indicates
  // that the reference minimization did not fully converge.
  if (!scan.points.empty()) {
    const auto min_it = std::min_element(
        scan.points.begin(), scan.points.end(),
        [](const UpperLimitPoint& a, const UpperLimitPoint& b) {
          return a.chi2 < b.chi2;
        });

    const double min_delta = min_it->chi2 - scan.best_chi2;
    if (min_delta < -1.0e-3) {
      std::cout
          << "\n[warning] A constrained scan point has chi2 below the "
          << "unconstrained reference by "
          << std::fixed << std::setprecision(4)
          << -min_delta
          << ". Check fit status and minimizer convergence.\n";
    }
  }

  std::cout
      << "\nUse the profile crossing together with residual inspection and "
      << "background/fit-window variations.\n\n";
}

} // namespace bwfit
