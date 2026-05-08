#pragma once

#include "FitConfig.h"

#include <memory>
#include <string>

#include <TH1F.h>

namespace bwfit {

double CalibrateEnergy(double branch_value, const FitConfig& cfg);
std::unique_ptr<TH1F> LoadSpectrum(const FitConfig& cfg);

} // namespace bwfit
