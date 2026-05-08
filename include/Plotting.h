#pragma once

#include "FitConfig.h"

#include <vector>

#include <TF1.h>
#include <TH1F.h>

namespace bwfit {

void DrawAndSaveFit(TH1F& hist, TF1& total_fit, TF1& no_interference_fit,
                    const FitConfig& cfg, const std::vector<double>& fitted_params,
                    double chi2, double ndf);

} // namespace bwfit
