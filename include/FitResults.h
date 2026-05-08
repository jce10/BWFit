// include/FitResults.h
#ifndef FITRESULTS_H
#define FITRESULTS_H

#include "FitConfig.h"
#include <TF1.h>

namespace bwfit {
void PrintFitParameters(const TF1& fit, const FitConfig& cfg);
}

#endif