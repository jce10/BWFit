#pragma once

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <TH1.h>
#include <TMatrixDSym.h>

namespace bwfit {

// All raw-fit positions and widths are in mm. Built-in widths are FWHM.
struct CalibrationParameter {
  double value = 0, min = 0, max = 0; // min == max fixes the parameter
};
struct CalibrationPeak {
  std::string name, shape = "gaussian";
  double energy_MeV = 0, energy_error_MeV = 0;
  bool anchor = true;
  CalibrationParameter yield, centroid, width;
};
struct CalibrationSpec {
  std::string rootfile, tree = "singles_tree", branch = "singles_fppos";
  std::string histogram; // if set, load TH1 instead of a tree
  std::string output_prefix = "outputs/calibration";
  double xmin_mm = -250, xmax_mm = 250;
  int bins = 2000, order = 1, background_order = 1;
  bool poisson = true;
  std::vector<CalibrationParameter> background; // Bernstein density, counts/mm
  std::vector<CalibrationPeak> peaks;
  int interference_first = -1, interference_second = -1;
  CalibrationParameter phase; // radians; raw-position orientation
};
struct CalibrationPeakResult {
  CalibrationPeak input;
  double yield = 0, yield_error = 0;
  double centroid_mm = 0, centroid_error_mm = 0;
  double fwhm_mm = 0, fwhm_error_mm = 0;
  double residual_keV = 0;
};
struct CalibrationResult {
  int order = 1, spectrum_status = -1, covariance_status = -1;
  int spectrum_ndf = 0, calibration_ndf = 0;
  double spectrum_chi2 = 0, spectrum_deviance = 0, calibration_chi2 = 0;
  // E(x) = a[0] + a[1]*x + a[2]*x*x, x in mm, E in MeV.
  std::array<double, 3> a{{0, 0, 0}};
  TMatrixDSym coefficient_covariance{3};
  TMatrixDSym spectrum_covariance;
  std::vector<double> spectrum_parameters, spectrum_errors;
  std::vector<CalibrationPeakResult> peaks;
  std::vector<std::string> warnings;
  double Energy(double x_mm) const;
  double EnergyError(double x_mm) const; // coefficient covariance only
  // First-order effective-variance GLS; preserves inter-peak correlations.
  // Detector systematics are not included or inferred from residuals.
};

// Custom shapes: positive unit-area density in 1/mm, (x, center, width).
// Register before ReadCalibrationSpec/RunCalibration. Additional shape
// parameters can be captured as fixed values in the callable.
using CalibrationShape = std::function<double(double, double, double)>;
void RegisterCalibrationShape(const std::string& name, CalibrationShape density);
CalibrationSpec ReadCalibrationSpec(const std::string& path);
std::unique_ptr<TH1> LoadRawSpectrum(const CalibrationSpec& spec);
CalibrationResult FitRawCalibrationSpectrum(TH1& histogram, const CalibrationSpec& spec);
CalibrationResult RunCalibration(const CalibrationSpec& spec);
void WriteCalibrationReports(TH1& histogram, const CalibrationSpec& spec,
                             const CalibrationResult& result);

} // namespace bwfit
