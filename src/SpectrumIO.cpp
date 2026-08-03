#include "SpectrumIO.h"

#include <iostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

namespace bwfit {

double CalibrateEnergy(double branch_value, const FitConfig& cfg) {
  return cfg.calib_slope * branch_value + cfg.calib_intercept;
}

int CalculateHistogramBins(const FitConfig& cfg) {
  const double fit_range = cfg.fitmax - cfg.fitmin;
  const double full_range = cfg.fullmax - cfg.fullmin;

  if (fit_range <= 0.0) {
    throw std::runtime_error(
        "CalculateHistogramBins: FIT_MAX must be greater than FIT_MIN.");
  }
  if (full_range <= 0.0) {
    throw std::runtime_error(
        "CalculateHistogramBins: FULL_MAX must be greater than FULL_MIN.");
  }
  if (cfg.fit_bins <= 0) {
    throw std::runtime_error(
        "CalculateHistogramBins: FIT_BINS must be greater than zero.");
  }

  const double bins_per_MeV =
      static_cast<double>(cfg.fit_bins) / fit_range;

  // Preserve the fitter's existing bin-count behavior exactly.
  const int nbins = static_cast<int>(full_range * bins_per_MeV);

  if (nbins <= 0) {
    throw std::runtime_error(
        "CalculateHistogramBins: calculated histogram bin count is not positive.");
  }

  return nbins;
}

double CalculateHistogramBinWidth(const FitConfig& cfg) {
  const int nbins = CalculateHistogramBins(cfg);
  return (cfg.fullmax - cfg.fullmin) / static_cast<double>(nbins);
}

std::unique_ptr<TH1F> LoadSpectrum(const FitConfig& cfg) {
  TFile infile(cfg.rootfile.c_str(), "READ");
  if (infile.IsZombie()) throw std::runtime_error("Could not open ROOT file: " + cfg.rootfile);

  TTree* tree = dynamic_cast<TTree*>(infile.Get(cfg.tree_name.c_str()));
  if (!tree) throw std::runtime_error("Could not find tree: " + cfg.tree_name);

  double branch_value = 0.0;
  tree->SetBranchAddress(cfg.branch_name.c_str(), &branch_value);

  const int nbins = CalculateHistogramBins(cfg);
  const double dE = CalculateHistogramBinWidth(cfg);

  auto hist = std::make_unique<TH1F>("fullhist", cfg.title.c_str(), nbins, cfg.fullmin, cfg.fullmax);
  hist->GetXaxis()->SetTitle("Excitation Energy (MeV)");
  hist->GetYaxis()->SetTitle(Form("Counts per %.0f keV", 1000.0 * dE));

  const Long64_t entries = tree->GetEntries();
  for (Long64_t i = 0; i < entries; ++i) {
    tree->GetEntry(i);
    hist->Fill(CalibrateEnergy(branch_value, cfg));
  }

  hist->SetDirectory(nullptr);
  return hist;
}

} // namespace bwfit
