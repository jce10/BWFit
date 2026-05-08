#include "FitConfig.h"
#include "FitModel.h"
#include "FitResults.h"
#include "SpectrumIO.h"
#include "Plotting.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TLine.h>
#include <TROOT.h>

using namespace bwfit;

namespace {

bool WantsLiveView(int argc, char* argv[]) {
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--view" || arg == "-v") return true;
  }
  return false;
}

void DrawLiveFitView(TH1& hist, TF1& total_fit, TF1& no_interference_fit, const FitConfig& cfg) {
  auto* live_canvas = new TCanvas("live_bw_fit", "Live BW Fit Viewer", 1400, 900);
  live_canvas->cd();

  hist.SetStats(0);
  hist.GetXaxis()->SetRangeUser(cfg.fullmin, cfg.fullmax);
  hist.GetXaxis()->SetTitle("Excitation Energy (MeV)");
  hist.GetYaxis()->SetTitle("Counts");
  hist.Draw("E");

  no_interference_fit.SetLineColor(kGray + 2);
  no_interference_fit.SetLineStyle(2);
  no_interference_fit.SetLineWidth(2);
  no_interference_fit.Draw("same");

  total_fit.SetLineColor(kBlue);
  total_fit.SetLineWidth(3);
  total_fit.Draw("same");

  const double ymin = hist.GetMinimum();
  const double ymax = hist.GetMaximum();

  TLine* fitmin_line = new TLine(cfg.fitmin, ymin, cfg.fitmin, ymax);
  fitmin_line->SetLineColor(kGreen + 2);
  fitmin_line->SetLineStyle(2);
  fitmin_line->SetLineWidth(2);
  fitmin_line->Draw("same");

  TLine* fitmax_line = new TLine(cfg.fitmax, ymin, cfg.fitmax, ymax);
  fitmax_line->SetLineColor(kGreen + 2);
  fitmax_line->SetLineStyle(2);
  fitmax_line->SetLineWidth(2);
  fitmax_line->Draw("same");

  live_canvas->Modified();
  live_canvas->Update();
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 2 || argc > 3) {
    std::cout << "Usage: ./bin/bwsingf config/angle_20.txt [--view]\n";
    return 1;
  }

  const bool view = WantsLiveView(argc, argv);

  int root_argc = 2;
  char web_off_arg[] = "--web=off";
  char* root_argv[] = {argv[0], web_off_arg, nullptr};
  std::unique_ptr<TApplication> app;
  if (view) {
    // Force ROOT to use the classic/native canvas instead of launching
    // the browser-based web canvas backend. This avoids Firefox profile
    // errors on systems where ROOT defaults to web display.
    app = std::make_unique<TApplication>("bwfit_app", &root_argc, root_argv);
    gROOT->SetBatch(kFALSE);
    gROOT->SetWebDisplay("off");
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

    TF1 full_fit("bw_fit_full", model, cfg.fullmin, cfg.fullmax, cfg.NumParams());
    full_fit.SetParameters(fitted.data());
    full_fit.SetNpx(cfg.graph_points);

    TF1 full_no_interf("bw_fit_full_nointerf", no_interference_model, cfg.fullmin, cfg.fullmax, cfg.NumParams());
    full_no_interf.SetParameters(fitted.data());
    full_no_interf.SetNpx(cfg.graph_points);

    double chi2 = fit.GetChisquare();
    double ndf = fit.GetNDF();

    DrawAndSaveFit(*hist, full_fit, full_no_interf, cfg, fitted, chi2, ndf);

    if (view) {
      DrawLiveFitView(*hist, full_fit, full_no_interf, cfg);
      std::cout << "\nLive canvas opened. Close the ROOT canvas/window to finish.\n";
      app->Run();
    }

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 2;
  }
}
