#include "Plotting.h"
#include "LineShapes.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
// #include <TColor.h>
// #include <TStyle.h>
#include <TString.h>

namespace bwfit {
namespace {

bool IsValidStateIndex(int index, const FitConfig& cfg) {
  return index >= 0 && index < static_cast<int>(cfg.states.size());
}

Color_t StateComponentColor(int index, const FitConfig& cfg) {
  if (index == cfg.trapped_index) return kGreen + 2;
  if (index == cfg.superrad_index) return kViolet + 1;
  return kRed;
}

int StateComponentLineWidth(int index, const FitConfig& cfg) {
  return (index == cfg.trapped_index || index == cfg.superrad_index) ? 3 : 2;
}

} // namespace

void DrawAndSaveFit(TH1F& hist, TF1& total_fit, TF1& no_interference_fit,
                    const FitConfig& cfg, const std::vector<double>& p,
                    double chi2, double ndf) {
  TCanvas canv("canv", "canv", 1400, 800);

  // --------------------------------------------------------------------------
  // Publication-style canvas layout:
  //   top pad    = spectrum + full fit + components
  //   bottom pad = residuals = data - fit
  //
  // The residual panel keeps the SAME full x-axis range as the spectrum panel,
  // but only draws residual points inside [cfg.fitmin, cfg.fitmax]. This keeps
  // the residuals visually aligned with the spectrum while matching the region
  // used by ROOT's "R" fit option.
  // --------------------------------------------------------------------------
  canv.cd();

  TPad* pad_spec = new TPad("pad_spec", "Spectrum",  0.0, 0.30, 1.0, 1.00);
  TPad* pad_res  = new TPad("pad_res",  "Residuals", 0.0, 0.00, 1.0, 0.30);

  pad_spec->SetBottomMargin(0.02);
  pad_spec->SetLeftMargin(0.12);
  pad_spec->SetRightMargin(0.04);

  pad_res->SetTopMargin(0.03);
  pad_res->SetBottomMargin(0.30);
  pad_res->SetLeftMargin(0.12);
  pad_res->SetRightMargin(0.04);

  pad_spec->Draw();
  pad_res->Draw();

  // -----------------------------
  // Top pad: spectrum and fit
  // -----------------------------
  pad_spec->cd();

  hist.SetLineColor(kBlack);
  hist.SetStats(0);

  // Let the residual panel carry the x-axis labels/title.
  hist.GetXaxis()->SetLabelSize(0.0);
  hist.GetXaxis()->SetTitleSize(0.0);
  hist.GetYaxis()->SetTitleSize(0.055);
  hist.GetYaxis()->SetLabelSize(0.045);
  hist.GetYaxis()->SetTitleOffset(0.85);

  hist.Draw();

  no_interference_fit.SetLineColor(kGray + 2);
  no_interference_fit.SetLineStyle(2);
  no_interference_fit.SetLineWidth(2);
  no_interference_fit.Draw("same");

  total_fit.SetLineColor(kBlue);
  total_fit.SetLineWidth(3);

  const int graph_points = cfg.graph_points;
  std::vector<double> E(graph_points), bg(graph_points);
  for (int i = 0; i < graph_points; ++i) {
    E[i] = cfg.fullmin + (cfg.fullmax - cfg.fullmin) * static_cast<double>(i) / (graph_points - 1);
    if (cfg.background_type == BackgroundType::Quadratic) {
      bg[i] = QuadraticBackground(E[i], p[cfg.BgEIndex()], p[cfg.BgA0Index()],
                                  p[cfg.BgA1Index()], p[cfg.BgA2Index()]);
    } else {
      bg[i] = LinearBackground(E[i], p[cfg.BgEIndex()], p[cfg.BgA0Index()], p[cfg.BgA1Index()]);
    }
  }
  TGraph bg_gr(graph_points, E.data(), bg.data());
  bg_gr.SetName("background");
  bg_gr.SetLineColor(kOrange + 7);
  bg_gr.SetLineStyle(3);
  bg_gr.SetLineWidth(2);
  bg_gr.Draw("L same");

  std::vector<std::unique_ptr<TF1>> components;
  for (int i = 0; i < static_cast<int>(cfg.states.size()); ++i) {
    TF1* comp = nullptr;
    if (cfg.states[i].is_bw && cfg.use_relativistic)
      comp = new TF1(Form("state_%d", i), BreitWignerRelRoot, cfg.fullmin, cfg.fullmax, 3);
    else if (cfg.states[i].is_bw)
      comp = new TF1(Form("state_%d", i), BreitWignerRoot, cfg.fullmin, cfg.fullmax, 3);
    else
      comp = new TF1(Form("state_%d", i), GaussianRoot, cfg.fullmin, cfg.fullmax, 3);

    comp->SetParameters(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
    comp->SetLineColor(StateComponentColor(i, cfg));
    comp->SetLineWidth(StateComponentLineWidth(i, cfg));
    comp->SetNpx(cfg.graph_points);
    comp->Draw("same");
    components.emplace_back(comp);
  }

  total_fit.Draw("same");

  // Neutral fit-range boundaries so they do not compete with the highlighted
  // trapped-state curve.
  const double ymax = hist.GetMaximum() * 1.1;
  TLine left_bound_spec(cfg.fitmin, 0.0, cfg.fitmin, ymax);
  TLine right_bound_spec(cfg.fitmax, 0.0, cfg.fitmax, ymax);
  left_bound_spec.SetLineColor(kGray + 1);
  right_bound_spec.SetLineColor(kGray + 1);
  left_bound_spec.SetLineStyle(7);
  right_bound_spec.SetLineStyle(7);
  left_bound_spec.Draw("same");
  right_bound_spec.Draw("same");

  TLegend* leg = pad_spec->BuildLegend();
  leg->Clear();
  const int ti = cfg.trapped_index;
  const int si = cfg.superrad_index;
  if (IsValidStateIndex(ti, cfg)) {
    leg->AddEntry(components[ti].get(),
                  Form("Trapped: E = %.3f MeV | #Gamma = %.3f MeV | Yield = %.0f",
                       p[3 * ti + 1], p[3 * ti + 2], p[3 * ti]),
                  "l");
  }
  if (IsValidStateIndex(si, cfg)) {
    leg->AddEntry(components[si].get(),
                  Form("Superrad: E = %.3f MeV | #Gamma = %.3f MeV | Yield = %.0f",
                       p[3 * si + 1], p[3 * si + 2], p[3 * si]),
                  "l");
  }

  const double chi2_ndf = (ndf > 0) ? chi2 / static_cast<double>(ndf) : 0.0;
  leg->AddEntry(&total_fit, Form("#chi^{2}/N = %.3f", chi2_ndf), "l");
  leg->Draw();

  // -----------------------------
  // Bottom pad: residuals
  // -----------------------------
  pad_res->cd();

  TGraphErrors residual_gr;
  residual_gr.SetName("residual_gr");
  residual_gr.SetTitle("");

  int ip = 0;
  double max_abs_residual = 0.0;

  for (int ibin = 1; ibin <= hist.GetNbinsX(); ++ibin) {
    const double x = hist.GetBinCenter(ibin);

    // Only draw residual points in the fitted region.
    if (x < cfg.fitmin || x > cfg.fitmax) continue;

    const double data = hist.GetBinContent(ibin);
    double err = hist.GetBinError(ibin);

    // If Sumw2 was not set, ROOT usually gives sqrt(N), but keep this safe.
    if (err <= 0.0) err = (data > 0.0) ? std::sqrt(data) : 1.0;

    const double fit = total_fit.Eval(x);
    const double res = data - fit;

    residual_gr.SetPoint(ip, x, res);
    residual_gr.SetPointError(ip, 0.0, err);

    if (std::fabs(res) + err > max_abs_residual) {
      max_abs_residual = std::fabs(res) + err;
    }

    ++ip;
  }

  if (max_abs_residual <= 0.0) max_abs_residual = 1.0;

  // Keep the same full x-axis range as the spectrum panel.
  const int residual_bins = hist.GetNbinsX();
  TH1F residual_frame("residual_frame", "", residual_bins, cfg.fullmin, cfg.fullmax);
  residual_frame.SetStats(0);
  residual_frame.GetXaxis()->SetTitle(hist.GetXaxis()->GetTitle());
  residual_frame.GetYaxis()->SetTitle("Residuals");

  residual_frame.GetXaxis()->SetTitleSize(0.12);
  residual_frame.GetXaxis()->SetLabelSize(0.10);
  residual_frame.GetYaxis()->SetTitleSize(0.10);
  residual_frame.GetYaxis()->SetLabelSize(0.08);
  residual_frame.GetYaxis()->SetTitleOffset(0.55);
  residual_frame.GetYaxis()->SetNdivisions(505);

  residual_frame.SetMinimum(-1.20 * max_abs_residual);
  residual_frame.SetMaximum( 1.20 * max_abs_residual);
  residual_frame.Draw();

  // Zero line only across the fitted region, where residuals are plotted.
  TLine zero_line(cfg.fitmin, 0.0, cfg.fitmax, 0.0);
  zero_line.SetLineColor(kBlue);
  zero_line.SetLineStyle(2);
  zero_line.Draw("same");

  // Matching neutral fit-range boundaries in the residual panel.
  TLine left_bound_res(cfg.fitmin, residual_frame.GetMinimum(), cfg.fitmin, residual_frame.GetMaximum());
  TLine right_bound_res(cfg.fitmax, residual_frame.GetMinimum(), cfg.fitmax, residual_frame.GetMaximum());
  left_bound_res.SetLineColor(kGray + 1);
  right_bound_res.SetLineColor(kGray + 1);
  left_bound_res.SetLineStyle(7);
  right_bound_res.SetLineStyle(7);
  left_bound_res.Draw("same");
  right_bound_res.Draw("same");

  residual_gr.SetMarkerStyle(20);
  residual_gr.SetMarkerSize(0.7);
  residual_gr.SetLineColor(kBlack);
  residual_gr.SetMarkerColor(kBlack);
  residual_gr.Draw("P same");

  canv.cd();

  TFile outfile(cfg.output_rootfile.c_str(), "RECREATE");
  canv.Write();
  hist.Write();
  total_fit.Write("total_fit");
  no_interference_fit.Write("no_interference_fit");
  residual_gr.Write("residual_gr");
  outfile.Close();
}

} // namespace bwfit
