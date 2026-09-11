#include "FitConfig.h"
#include "FitModel.h"
#include "FitResults.h"
#include "SpectrumIO.h"
#include "Plotting.h"
#include "UpperLimit.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

using namespace bwfit;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "Usage: ./bin/bwsingf config/angle_20.txt\n";
    return 1;
  }

  try {
    FitConfig cfg = ReadConfig(argv[1]);
    PrintConfig(cfg);

    auto hist = LoadSpectrum(cfg);

    FitModel model(cfg, true);
    FitModel no_interference_model(cfg, false);

    std::vector<double> params = InitialParameterArray(cfg);
    TF1 fit("bw_fit", model, cfg.fitmin, cfg.fitmax, cfg.NumParams());
    fit.SetParameters(params.data());
    fit.SetNumberFitPoints(cfg.graph_points);
    fit.SetNpx(cfg.graph_points);
    ConfigureParameterLimits(fit, cfg);

    std::cout << "Fitting...\n";
    TFitResultPtr result = hist->Fit(&fit, "RS", "", cfg.fitmin, cfg.fitmax);
    const int fit_status = static_cast<int>(result);
    std::cout << "Fit status = " << fit_status << "\n";
    std::vector<double> fitted(cfg.NumParams());
    for (int i = 0; i < cfg.NumParams(); ++i) {
      fitted[i] = fit.GetParameter(i);
    }

    std::cout << "\nChi2/NDF = "
              << fit.GetChisquare() / fit.GetNDF()
              << "\n\n";

    PrintFitParameters(fit, cfg);
    PrintFitDiagnostics(*result, cfg);

    char save_json = 'n';
    std::cout << "Save full fit summary to JSON? (y/n): ";
    std::cin >> save_json;
    if (save_json == 'y' || save_json == 'Y') {
      std::string json_name;
      std::cout << "Enter JSON filename: ";
      std::cin >> json_name;

      WriteFitSummaryJSON(
          fit,
          *result,
          cfg,
          json_name,
          fit_status,
          argv[1]);
    }
    // ---------------------------------------------------------------------- //
    // Profile-likelihood upper-limit scan for the superradiant normalization.
    //
    // The scan variable is the fitter's normalization parameter in
    // counts*MeV. UpperLimit.cpp also prints the corresponding integral counts
    // using integral = normalization / histogram bin width.
    //
    // Current grid:
    //   0.0 to 30.0 counts*MeV in steps of 0.1 counts*MeV
    //
    // With approximately 12-keV bins, each step is roughly eight counts.
    // Increase the upper endpoint if Delta chi2 = 2.71 is not reached.
    // ---------------------------------------------------------------------- //
    // const double scan_yield_min = 0.0;
    // const double scan_yield_max = 30.0;
    // const int scan_points = 301;

    // const std::vector<double> yields =
    //     MakeYieldGrid(scan_yield_min, scan_yield_max, scan_points);

    // std::cout << "\nRunning superradiant upper-limit scan from "
    //           << scan_yield_min << " to " << scan_yield_max
    //           << " counts*MeV using " << scan_points << " points...\n";
    // // Match the scan model to the nominal fit configuration.
    // const bool scan_interference = cfg.use_interference;

    // const UpperLimitScanResult scan =
    //     ScanUpperLimit(
    //         *hist,
    //         cfg,
    //         fitted,
    //         yields,
    //         scan_interference);

    // PrintUpperLimitScan(scan);

    // ---------------------------------------------------------------------- //
    TF1 full_fit(
        "bw_fit_full",
        model,
        cfg.fullmin,
        cfg.fullmax,
        cfg.NumParams());
    full_fit.SetParameters(fitted.data());
    full_fit.SetNpx(cfg.graph_points);

    TF1 full_no_interf(
        "bw_fit_full_nointerf",
        no_interference_model,
        cfg.fullmin,
        cfg.fullmax,
        cfg.NumParams());
    full_no_interf.SetParameters(fitted.data());
    full_no_interf.SetNpx(cfg.graph_points);
    const double chi2 = fit.GetChisquare();
    const double ndf = fit.GetNDF();

    DrawAndSaveFit(
        *hist,
        full_fit,
        full_no_interf,
        cfg,
        fitted,
        chi2,
        ndf);

    return 0;

  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 2;
  }
}
