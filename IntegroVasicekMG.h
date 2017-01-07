#ifndef __INTEGROVASICEKMG_H_INCLUDED__
#define __INTEGROVASICEKMG_H_INCLUDED__

#include <cmath>
#include <vector>
#include <complex>
#include "FunctionalUtilities.h"

namespace creditriskutilities {
    struct loan{
        double pd;
        //std::function<Complex(const Complex&)> lgdCF;//characteristic function
        double exposure;
        std::vector<double> w;
        loan(double pd_, double exposure_, std::vector<double>&& w_){
            pd=pd_;
            exposure=exposure_;
            w=w_;
        }
    };
    /** u is typically complex*/
    template<typename Number, typename Loans, typename LGDCF, typename Lambda, typename Q>
    auto logLPMCF(const Number &u, const std::vector<Loans>& loans, const LGDCF& lgdCF, const Lambda& lambda, const Q& q){
        int n=loans.size();
        int m=loans[0].w.size();
        auto upperU=-lambda*u; //liquidity risk..u*lambda*i
        upperU=(exp(upperU)-1.0)*q;//.multiply(new Complex(0, 1));//q*(exp(i*u*lambda)-1)
        upperU=upperU-u; //liquidity..u*i+q*(exp(u*lambda*i)-1)
        return futilities::for_each_parallel(0, m, [&](const auto& indexM){
            futilities::sum(0, n, [&](const auto& index){
                return (lgdCF(upperU, loans[index])-1.0)*loans[index].pd*loans[index].w[indexM];
            });
        });
    }


    template<typename Number, typename LoanSize, typename Lambda, typename Theta, typename Sigma, typename T, typename X0>
    auto lgdCF(const std::complex<Number> &u, const LoanSize &l, const Lambda &lambda,const Theta &theta, const Sigma &sigma, const T &t, const X0 &x0){/*I think this is a CIR characteristic function*/
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
        //assert(y0.size()==alpha.size());

        return futilities::for_each_parallel_copy(y0, [&](const auto& val, const auto& index){
            return (val-1)*helpComputeMoments(alpha[index], tau);
        });
    }
    /**
    Computes the variance/covariance of the integral of a multivariate vasicek process
    */
    template<typename Alpha, typename Tau, typename Sigma, typename Rho>
    auto computeVarianceVasicek(const std::vector<Alpha>& alpha, const std::vector<Sigma>& sigma, const std::vector<std::vector<Rho> >& rho,  const Tau& tau){
        //assert(alpha.size()==sigma.size()&&rho.size()==alpha.size()&&rho[0].size()==alpha.size());
        int rowLength=alpha.size();
        return futilities::for_each_parallel(0, rowLength, [&](const auto& indexI){
            auto ai=helpComputeMoments(alpha[indexI], tau);
            return futilities::for_each_parallel(0, rowLength, [&](const auto& indexJ){
                auto aj=helpComputeMoments(alpha[indexJ], tau);
                return crossMultiply(rho[indexI][indexJ], sigma[indexI], sigma[indexJ], alpha[indexI], alpha[indexJ])*(tau-ai-aj+helpComputeMoments(alpha[indexI]+alpha[indexJ], tau));
            });
        });
    }

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
