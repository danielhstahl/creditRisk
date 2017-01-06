#ifndef __INTEGROVASICEKMG_H_INCLUDED__
#define __INTEGROVASICEKMG_H_INCLUDED__

#include <cmath>
#include <vector>
#include <complex>
#include "FunctionalUtilities.h"

namespace CreditRiskUtilities {
    template<typename Alpha, typename Tau>
    auto helpComputeMoments(const Alpha& alpha, const Tau& tau){ //hleper function since called so much
        return (1-exp(-alpha*tau))/alpha;
    }
    template<typename Rho, typename Sigma, typename Alpha>
    auto crossMultiply(const Rho& rho, const Sigma& sigma1, const Sigma& sigma2, const Alpha& alpha1, const Alpha& alpha2){
        return (rho*sigma1*sigma2)/(alpha1*alpha2);
    }

    /**
    Computes the expectation of the integral of a multivariate vasicek process
    */
    template<typename Number, typename Alpha, typename Tau>
    auto computeExpectationVasicek(const std::vector<Number>& y0, const std::vector<Alpha>& alpha, const Tau& tau){
        assert(y0.size()==alpha.size());
        return futilities::for_each_parallel(std::move(alpha), [&](val, index){
            return (y0[index]-1)*helpComputeMoments(val, tau);
        });
    }
    /**
    Computes the variance/covariance of the integral of a multivariate vasicek process
    */
    template<typename Alpha, typename Tau, typename Sigma, typename Rho>
    auto computeVarianceVasicek(const std::vector<Alpha>& alpha, const std::vector<Sigma>& sigma, const std::vector<std::vector<Rho> >& rho,  const Tau& tau){
        assert(alpha.size()==sigma.size()&&rho.size()==rho[0].size()&&rho.size()==alpha.size());
        return futilities::for_each_parallel(std::move(rho), [&](const auto& row, const auto& indexI){
            auto ai=helpComputeMoments(alpha[indexI], tau);
            return futilities::for_each_parallel(std::move(row), [&](const auto& val, const auto& indexJ){
                auto aj=helpComputeMoments(alpha[indexJ], tau);
                return crossMultiply(val, sigma[indexI], sigma[indexJ], alpha[indexI], alpha[indexJ])*(tau-ai-aj+helpComputeMoments(alpha[indexI]+alpha[indexJ], tau));
            });
        });
    }


    /*template<typename number>
    void computeVasicekMGF(std::vector<number> &expectation, std::vector< std::vector<number> > &variance, const std::vector<double> &alpha, const std::vector<double> &sigma, const std::vector<std::vector<double> > &rho, const std::vector<double> &y0, double tau){
        int m=alpha.size();
        for(int i=0; i<m; i++){
            double ai=helpComputeMoments(alpha[i], tau);
            expectation.push_back((y0[i]-1)*ai+tau);
            std::vector<double> row=std::vector<double>(m);
            for(int j=0; j<m; j++){
                double aj=helpComputeMoments(alpha[j], tau);
                double helpVarij=(rho[i][j]*sigma[i]*sigma[j]/(alpha[i]*alpha[j]))*(tau-ai-aj+(1-std::exp(-(alpha[i]+alpha[j])*tau))/(alpha[i]+alpha[j])); //difra page 10
                row[j]=helpVarij;
            }
            variance.push_back(row);
        }
    }*/
    template<typename Expectation, typename Variance>
    auto executeVasicekMGF(const std::vector<std::complex> &v, const std::vector<Expectation>& expectation , const std::vector< std::vector<Variance> >& variance){

        auto m=expectation.size();
        return exp(futilities::sum(expectation, [&](const auto& val, const auto& index){
            return v[index]*val;
        })+futilities::sum(0, m, [&](const auto& indexI){
            futilities::sum(variance[indexI], [&](const auto& val, const auto& indexJ){
                return v[indexI]*v[indexJ]*val;
            });
        })*.5);
       /* for(int i=0; i<m; i++){
            el=el.add(v[i].multiply(expectation[i]));
            for(int j=0; j<m; j++){
                var=var.add(v[i].multiply(v[j]).multiply(variance[i][j]));
            }
        }
        var=var.multiply(.5).add(el);
        Complex phi=var.exp();
        return phi;*/
    }


}

#endif
