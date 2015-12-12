#ifndef __INTEGROVASICEKMG_H_INCLUDED__
#define __INTEGROVASICEKMG_H_INCLUDED__

#include <cmath>
#include <vector>
#include "Complex.h"

auto helpComputeMoments(
    const auto&, 
    const auto&
);
template<typename number>
void computeVasicekMGF( //computes mu and Sigma for multidimensional integral of vasicek
    std::vector<number>&, //expectation (begins with empty vector)
    std::vector< std::vector<number> >&, //variance (begins with empty vector)
    const std::vector<double>&,//alpha
    const std::vector<double>&, //sigma
    const std::vector<std::vector<double> >&, //rho
    const std::vector<double>&, //y0
    double //tau
);
template<typename number>
Complex executeVasicekMGF(
    const std::vector<Complex>&, //v in E[e^vX]
    const std::vector<number>&, //expectation
    const std::vector< std::vector<number> >& //variance
);
#include "IntegroVasicekMG.hpp"
#endif
