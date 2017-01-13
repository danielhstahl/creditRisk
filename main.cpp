#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOost.h"
#include <cmath>
#include "CreditUtilities.h"
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson
#include "FunctionalUtilities.h"

template<typename Container, typename Range>
void printJson(const Container& myContainer, const Range& mn, const Range& dx){
	std::cout<<"{\"y\":[";
	for(const auto& val:myContainer){
		std::cout<<val<<",";
	}
	std::cout<<"],\"xmin\":"<<mn<<",\"dx\":"<<dx<<"}";
}
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
		loan(){
    }
};

template<typename RAND>
auto getP(double minP, double maxP,const RAND& rand){
	return maxP*rand()/RAND_MAX+minP;
}
template<typename RAND>
auto getExposure(double minLoan, double maxLoan, const RAND& rand){
	return (maxLoan-minLoan)*rand()/RAND_MAX+minLoan;
}
template<typename RAND>
auto getNormalizedWeights(int numWeights, const RAND& rand){
	double totalW=0;
	return futilities::for_each_parallel(futilities::for_each_parallel(0, numWeights, [&](const auto& index){
		auto w=1.0*rand();
		totalW+=w;
		return w;
	}), [&](const auto& val, const auto& index){
		return val/totalW;
	});
}
int main(int argc, char* argv[]){
	int xNum=1024;
	int uNum=256;
  int n=100000;//sensible default values
	int m=1;
	double tau=1;
	double qUnscaled=.05;
	double lambdaUnscaled=.05;
	const double minLoan=10000;
	const double maxLoan=50000;
	const double minP=.0001;
	double maxP=.09;
	std::vector<double> alpha(m, .2);
	std::vector<double> sigma(m, .3);
	std::vector<double> y0(m, .5);
	std::vector< std::vector<double> > rho(m, std::vector<double> (m, 1.0));
	//These are parameters for the LGD CF function
  const double alphL=.2;
	const double bL=.5;
	const double sigL=.2;
	srand(5);
	if(argc>1){
		rapidjson::Document parms;
		parms.Parse(argv[1]);//yield data
		if(parms.FindMember("xNum")!=parms.MemberEnd()){
			xNum=parms["xNum"].GetInt();
		}
		if(parms.FindMember("uNum")!=parms.MemberEnd()){
			uNum=parms["uNum"].GetInt();
		}
		if(parms.FindMember("n")!=parms.MemberEnd()){
			n=parms["n"].GetInt();
		}
		if(parms.FindMember("q")!=parms.MemberEnd()){
			qUnscaled=parms["q"].GetDouble();
		}
		if(parms.FindMember("lambda")!=parms.MemberEnd()){
			lambdaUnscaled=parms["lambda"].GetDouble();
		}
		if(parms.FindMember("alpha")!=parms.MemberEnd()){
			alpha[0]=parms["alpha"].GetDouble();
		}
		if(parms.FindMember("sigma")!=parms.MemberEnd()){
			sigma[0]=parms["sigma"].GetDouble();
		}
		if(parms.FindMember("t")!=parms.MemberEnd()){
			tau=parms["t"].GetDouble();
		}
		if(parms.FindMember("maxP")!=parms.MemberEnd()){
			maxP=parms["maxP"].GetDouble();
		}
		if(parms.FindMember("x0")!=parms.MemberEnd()){
			y0[0]=parms["x0"].GetDouble();
		}
	}
	const auto loans=futilities::for_each_parallel(0, n, [&](const auto& index){
		return loan(getP( minP, maxP, rand), getExposure(minLoan, maxLoan, rand), getNormalizedWeights(m, rand));
	});

	const double expectedTotalExposure=(minLoan+.5*(maxLoan-minLoan))*n;
	const double lambda=lambdaUnscaled*expectedTotalExposure; //proxy for n*exposure*lambda
	const double q=qUnscaled/lambda;
	const double xmax=0;
	const double xmin=-expectedTotalExposure*bL*maxP*.5*5;//5 is arbitrary

	const auto expectation=creditutilities::computeExpectationVasicek(y0, alpha, tau);
	const auto variance=creditutilities::computeVarianceVasicek(alpha, sigma, rho, tau);

	const auto density=fangoost::computeInv(xNum, uNum,  xmin, xmax, [&](const auto& u){
		return creditutilities::executeVasicekMGF(
				creditutilities::logLPMCF(
					creditutilities::getUpperU(u, lambda, q),
					loans,
					n, m, 
					[&](const auto& u, const auto& l){
						return creditutilities::lgdCF(u, l.exposure, alphL, bL, sigL, tau, bL);
					},
					[](const auto& loan){
						return loan.pd;
					},
					[](const auto& loan, const auto& index){
						return loan.w[index];
					}
				),
				expectation,
				variance
			);
	});
	const auto dx=fangoost::computeDX(xNum, xmin, xmax);
	printJson(density, xmin, dx);
	
}
