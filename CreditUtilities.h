#ifndef __CREDITUTILITIES_H_INCLUDED__
#define __CREDITUTILITIES_H_INCLUDED__

#include <cmath>
#include <vector>
#include <complex>
#include "FunctionalUtilities.h"

namespace creditutilities {
    

    template<typename Number, typename Lambda, typename Q>
    auto getUpperU(const Number& u, const Lambda& lambda, const Q& q){
        return -(exp(-u*lambda)-1.0)*q-u;
    }

    /** This function returns the exponent for the credit risk characteristic function: namely sum_j p_j Y (phi_j(u)-1).  Note that to do liquidity risk, pass "getUpperU"'s results to "u". */
    template<typename Number, typename Loans, typename Incr, typename GetW, typename GetPD, typename LGDCF>
    auto logLPMCF(const Number &u, const std::vector<Loans>& loans, const Incr& m, const LGDCF& lgdCF, const GetPD& getPD, const GetW& getW){
        return futilities::for_each_parallel(0, m, [&](const auto& indexM){
            Incr n=loans.size();
            return futilities::sum(0, n, [&](const auto& index){
                return (lgdCF(u, loans[index])-1.0)*getPD(loans[index])*getW(loans[index], indexM);
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

    template<typename Alpha, typename Tau>
    auto helpComputeMoments(const Alpha& alpha, const Tau& tau){ //hleper function since called so much
        return (1-exp(-alpha*tau))/alpha;
    }
    template<typename Rho, typename Sigma, typename Alpha>
    auto crossMultiply(const Rho& rho, const Sigma& sigma1, const Sigma& sigma2, const Alpha& alpha1, const Alpha& alpha2){
        return (rho*sigma1*sigma2)/(alpha1*alpha2);
    }

    /**
    Computes the expectation of the integral of a multivariate vasicek process with mean 1
    */
    template<typename Number, typename Alpha, typename Tau>
    auto computeExpectationVasicek(const std::vector<Number>& y0, const std::vector<Alpha>& alpha, const Tau& tau){
        return futilities::for_each_parallel_copy(y0, [&](const auto& val, const auto& index){
            return (val-1)*helpComputeMoments(alpha[index], tau);
        });
    }
    /**
    Computes the variance/covariance of the integral of a multivariate vasicek process
    */
    template<typename Alpha, typename Tau, typename Sigma, typename Rho>
    auto computeVarianceVasicek(const std::vector<Alpha>& alpha, const std::vector<Sigma>& sigma, const std::vector<std::vector<Rho> >& rho,  const Tau& tau){
        int rowLength=alpha.size();
        return futilities::for_each_parallel(0, rowLength, [&](const auto& indexI){
            auto ai=helpComputeMoments(alpha[indexI], tau);
            return futilities::for_each_parallel(0, rowLength, [&](const auto& indexJ){
                auto aj=helpComputeMoments(alpha[indexJ], tau);
                return crossMultiply(rho[indexI][indexJ], sigma[indexI], sigma[indexJ], alpha[indexI], alpha[indexJ])*(tau-ai-aj+helpComputeMoments(alpha[indexI]+alpha[indexJ], tau));
            });
        });
    }
    /**Computes the expectation of a the exponential of a weighted combination of the multidemensional integrated vasicek process*/
    template<typename Expectation, typename Variance, typename Number>
    auto executeVasicekMGF(const std::vector<std::complex<Number> > &v, const std::vector<Expectation>& expectation , const std::vector< std::vector<Variance> >& variance){

        int m=expectation.size();
        return exp(futilities::sum(0, m, [&](const auto& index){
            return v[index]*expectation[index];
        })+futilities::sum(0, m, [&](const auto& indexI){
            return futilities::sum(0, m, [&](const auto& indexJ){
                return v[indexI]*v[indexJ]*variance[indexI][indexJ];
            });
        })*.5);

    }


}

#endif
