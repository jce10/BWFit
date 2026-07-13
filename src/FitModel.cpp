#include "FitModel.h"
#include "LineShapes.h"

#include <cmath>
#include <iostream>

namespace bwfit {

FitModel::FitModel(const FitConfig& cfg, bool interference_enabled)
    : cfg_(cfg), interference_enabled_(interference_enabled) {
  num_states_ = static_cast<int>(cfg_.states.size());
  num_params_ = cfg_.NumParams();
}

double FitModel::operator()(double* var_arr, double* param_arr) {
  

  const double E = var_arr[0];
  double return_val = 0.0;

  std::vector<double> fac(num_states_);
  std::vector<double> M(num_states_);
  std::vector<double> G(num_states_);

  for (int i = 0; i < num_states_; ++i) {
    fac[i] = param_arr[3 * i + 0];
    M[i]   = param_arr[3 * i + 1];
    G[i]   = param_arr[3 * i + 2];

    if (cfg_.states[i].is_bw && cfg_.use_relativistic)
      return_val += BreitWignerRel(E, M[i], G[i], fac[i]);

    else if (cfg_.states[i].is_bw && !cfg_.use_relativistic)
      return_val += BreitWigner(E, M[i], G[i], fac[i]);

    else
      return_val += Gaussian(E, M[i], G[i], fac[i]);
  }


  const double phase = param_arr[cfg_.PhaseIndex()];
  const double bg_e  = param_arr[cfg_.BgEIndex()];
  const double bg_a0 = param_arr[cfg_.BgA0Index()];
  const double bg_a1 = param_arr[cfg_.BgA1Index()];


  if (cfg_.background_type == BackgroundType::Quadratic) {
    const double bg_a2 = param_arr[cfg_.BgA2Index()];
    return_val += QuadraticBackground(E, bg_e, bg_a0, bg_a1, bg_a2);
  } 
  else {
    return_val += LinearBackground(E, bg_e, bg_a0, bg_a1);
  }

  if (interference_enabled_ && cfg_.use_interference) {
    const int TI = cfg_.trapped_index;
    const int SI = cfg_.superrad_index;

    if (TI >= 0 && TI < num_states_ && SI >= 0 && SI < num_states_) {
      double X = 0.0, Y = 0.0, F = 0.0;
      
      if (!cfg_.use_relativistic) {
        X = (E - M[TI]) * (E - M[SI]) + G[TI] * G[SI] / 4.0;

        Y = (G[TI] / 2.0) * (E - M[SI]) - (G[SI] / 2.0) * (E - M[TI]);
        
        F = 2.0 * std::sqrt(fac[TI] * fac[SI] * G[TI] * G[SI]) / (2.0 * M_PI) *
            (X * std::cos(phase) + Y * std::sin(phase)) / (X * X + Y * Y);
      } 
      else {
        X = (E * E - M[TI] * M[TI]) * (E * E - M[SI] * M[SI]) +
            M[TI] * M[SI] * G[TI] * G[SI];

        Y = M[TI] * G[TI] * (E * E - M[SI] * M[SI]) -
            M[SI] * G[SI] * (E * E - M[TI] * M[TI]);

        const double kT = KFactor(M[TI], G[TI]);
        const double kS = KFactor(M[SI], G[SI]);

        F = 2.0 * std::sqrt(fac[TI] * fac[SI] * kT * kS) *
            (X * std::cos(phase) + Y * std::sin(phase)) / (X * X + Y * Y);
      }
      return_val += F;
    }
  }

  return return_val;
}

std::vector<double> InitialParameterArray(const FitConfig& cfg) {
  std::vector<double> p(cfg.NumParams());
  for (int i = 0; i < static_cast<int>(cfg.states.size()); ++i) {
    p[3 * i + 0] = cfg.states[i].fac;
    p[3 * i + 1] = cfg.states[i].M;
    p[3 * i + 2] = cfg.states[i].G;
  }

  p[cfg.PhaseIndex()] = cfg.phase.value;
  p[cfg.BgEIndex()]   = cfg.bg_e.value;
  p[cfg.BgA0Index()]  = cfg.bg_a0.value;
  p[cfg.BgA1Index()]  = cfg.bg_a1.value;
  if (cfg.background_type == BackgroundType::Quadratic) {
    p[cfg.BgA2Index()] = cfg.bg_a2.value;
  }
  return p;
}

void ConfigureParameterLimits(TF1& fit, const FitConfig& cfg) {
  for (int i = 0; i < 3 * static_cast<int>(cfg.states.size()); ++i) {
    const int state = i / 3;
    switch (i % 3) {
      case 0:
        fit.SetParLimits(i, 0.0, 1.0E8);
        break;
      case 1:
        fit.SetParLimits(i, cfg.fitmin, cfg.fitmax);
        if (!cfg.states[state].fit_M) fit.FixParameter(i, cfg.states[state].M);
        break;
      case 2: {
        const double min_G = cfg.states[state].G * (1.0 - cfg.width_perc);
        const double max_G = cfg.states[state].G * (1.0 + cfg.width_perc);
        fit.SetParLimits(i, min_G, max_G);
        if (!cfg.states[state].fit_G) fit.FixParameter(i, cfg.states[state].G);
        break;
      }
    }
  }

  if (cfg.fit_phase) {
    fit.SetParLimits(cfg.PhaseIndex(), cfg.phase.min, cfg.phase.max);
  } else {
    fit.FixParameter(cfg.PhaseIndex(), cfg.phase.value);
  }
  fit.SetParLimits(cfg.BgEIndex(), cfg.bg_e.min, cfg.bg_e.max);
  fit.SetParLimits(cfg.BgA0Index(), cfg.bg_a0.min, cfg.bg_a0.max);
  fit.SetParLimits(cfg.BgA1Index(), cfg.bg_a1.min, cfg.bg_a1.max);
  if (cfg.background_type == BackgroundType::Quadratic) {
    fit.SetParLimits(cfg.BgA2Index(), cfg.bg_a2.min, cfg.bg_a2.max);
  }
}

} // namespace bwfit
