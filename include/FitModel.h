#pragma once

#include "FitConfig.h"

#include <vector>

#include <TF1.h>

namespace bwfit {

class FitModel {
public:
  FitModel(const FitConfig& cfg, bool interference_enabled);
  double operator()(double* x, double* p);

  int NumParams() const { return num_params_; }

private:
  const FitConfig& cfg_;
  bool interference_enabled_ = true;
  int num_states_ = 0;
  int num_params_ = 0;
};

std::vector<double> InitialParameterArray(const FitConfig& cfg);
void ConfigureParameterLimits(TF1& fit, const FitConfig& cfg);

} // namespace bwfit
