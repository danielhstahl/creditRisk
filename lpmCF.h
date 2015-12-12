#ifndef __LPMCF_H_INCLUDED__
#define __LPMCF_H_INCLUDED__


#include <cmath>
#include "Complex.h"
#include "Asset.h"
#include <vector>

std::vector<Complex> logLPMCF(
    Complex&, //u
    std::vector<loan>&, //loans
    const auto&, //lgd function
    double, //lambda
    double //q
);
#include "lpmCF.hpp"
#endif
