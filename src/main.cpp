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
    std::cout << "Fit status = " << static_cast<int>(result) << "\n";

    std::vector<double> fitted(cfg.NumParams());
    for (int i = 0; i < cfg.NumParams(); ++i) fitted[i] = fit.GetParameter(i);

    std::cout << "\nChi2/NDF = " << fit.GetChisquare() / fit.GetNDF() << "\n\n";
    PrintFitParameters(fit, cfg);

    // -------------------------------------------------------------------------
// Upper limit scan for the superradiant state
// -------------------------------------------------------------------------

    // std::vector<double> yields = {
    //   0.0, 0.1, 0.5, 1.0, 2.0,
    //   5.0, 10.0, 20.0, 50.0, 100.0
    // };
    // auto scan = ScanUpperLimit(*hist, cfg, fitted, yields, false);
    // PrintUpperLimitScan(scan);

    TF1 full_fit("bw_fit_full", model, cfg.fullmin, cfg.fullmax, cfg.NumParams());
    full_fit.SetParameters(fitted.data());
    full_fit.SetNpx(cfg.graph_points);

    TF1 full_no_interf("bw_fit_full_nointerf", no_interference_model, cfg.fullmin, cfg.fullmax, cfg.NumParams());
    full_no_interf.SetParameters(fitted.data());
    full_no_interf.SetNpx(cfg.graph_points);

    double chi2 = fit.GetChisquare();
    double ndf = fit.GetNDF();

    DrawAndSaveFit(*hist, full_fit, full_no_interf, cfg, fitted, chi2, ndf);

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 2;
  }
}
