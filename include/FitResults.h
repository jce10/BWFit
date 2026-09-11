// include/FitResults.h
#ifndef FITRESULTS_H
#define FITRESULTS_H

#include "FitConfig.h"
#include <string>
#include <TF1.h>
#include <TFitResult.h>

namespace bwfit {

void PrintFitParameters(const TF1& fit, const FitConfig& cfg);

void PrintFitDiagnostics(
    const TFitResult& result,
    const FitConfig& cfg
);

void WriteFitSummaryJSON(
    const TF1& fit,
    const TFitResult& result,
    const FitConfig& cfg,
    const std::string& output_json,
    int fit_status,
    const std::string& config_path = ""
);

}

#endif
