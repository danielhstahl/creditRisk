#ifndef __CREDITUTILITIES_H_INCLUDED__
#define __CREDITUTILITIES_H_INCLUDED__

#include <cmath>
#include <vector>
#include <complex>
#include "FunctionalUtilities.h"


/**This is not in its own repo because it is choices that must be made per implementation of each credit model*/
namespace creditutilities {
    

    template<typename Number, typename Lambda, typename Q>
    auto getLiquidityRisk(const Number& u, const Lambda& lambda, const Q& q){
        return -(exp(-u*lambda)-1.0)*q-u;
    }


    
    /** This function returns the exponent for the credit risk characteristic function: namely sum_j p_j Y (phi_j(u)-1).  Note that to do liquidity risk, pass "getLiquidityRisk"'s results to "u". */
    template<typename Number, typename Loans, typename Incr, typename GetW, typename GetPD, typename LGDCF>
    auto logLPMCF(const Number &u, const std::vector<Loans>& loans, const Incr& m, const LGDCF& lgdCF, const GetPD& getPD, const GetW& getW){
        return futilities::for_each_parallel(0, m, [&](const auto& indexM){
            return futilities::sum(loans, [&](const auto& loan, const auto& index){
                return (lgdCF(u, loan)-1.0)*getPD(loan)*getW(loan, indexM);
            });
        });
    }



    /**Characteristic function for LGD.  Follows CIR process.  U is typically complex*/
    template<typename Number, typename LoanExposure, typename Lambda, typename Theta, typename Sigma, typename T, typename X0>
    auto lgdCF(const Number &u, const LoanExposure &l, const Lambda &lambda,const Theta &theta, const Sigma &sigma, const T &t, const X0 &x0){
        auto expt=exp(-lambda*t);
        auto sigL=-sigma*sigma/(2*lambda);
        auto uu=u*l;
        auto uP=uu*(1-expt)*sigL+1.0;
        return exp((uu*expt*x0)/uP)*pow(uP, theta/sigL);
    }
};
#endif
