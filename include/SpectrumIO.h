#pragma once

#include "FitConfig.h"

#include <memory>

#include <TH1F.h>

namespace bwfit {

double CalibrateEnergy(double branch_value, const FitConfig& cfg);

// Return the number of bins used to construct the full histogram.
int CalculateHistogramBins(const FitConfig& cfg);

// Return the exact width of the histogram bins in MeV/bin.
double CalculateHistogramBinWidth(const FitConfig& cfg);

std::unique_ptr<TH1F> LoadSpectrum(const FitConfig& cfg);

} // namespace bwfit
