#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOost.h"
#include <cmath>
#include "IntegroVasicekMG.h"
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson
#include "FunctionalUtilities.h"

template<typename Container, typename Range>
void printJson(const Container& myContainer, const Range& mn, const Range& dx){
	std::cout<<"{y:[";
	for(const auto& val:myContainer){
		std::cout<<val<<",";
	}
	std::cout<<"],xmin:"<<mn<<",dx:"<<dx<<"}\\n";
}

int main(int argc, char* argv[]){
	int xNum=1024;
	int uNum=256;
  int n=100000;//sensible defaults
	int m=1;
	double tau=1;
	std::vector<double> alpha=std::vector<double>(m);
	std::vector<double> sigma=std::vector<double>(m);
	std::vector<double> y0=std::vector<double>(m);
	std::vector< std::vector<double> > rho(m, std::vector<double> (m));
	

	//These are parameters for the LGD CF function
  const double alphL=.2;
	const double bL=.5;
	const double sigL=.2;
	srand(5);
	double maxP=.09;
  //double totalwP=0;
	double q=.05;
	double lambda=.05;
	const double minLoan=10000;
	const double maxLoan=50000;
	const double diffLoan=maxLoan-minLoan;
	const double minP=.0001;
	rapidjson::Document parms;
	parms.Parse(argv[0]);//yield data
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
		q=parms["q"].GetDouble();
	}
	if(parms.FindMember("lambda")!=parms.MemberEnd()){
		lambda=parms["lambda"].GetDouble();
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
	auto getNormalizedWeights=[&](int numWeights){
		double totalW=0;
		return futilities::for_each(futilities::for_each_parallel(0, numWeights, [&](const auto& index){
			auto w=1.0*rand();
			totalW+=w;
			return w;
		}), [&](const auto& val, const auto& index){
			return val/totalW;
		});
	};
	auto getP=[&](double maxP){
		return maxP*rand()/RAND_MAX+minP;
	};
	double totalExposure=0;
	auto loans=futilities::for_each_parallel(0, n, [&](const auto& index){
		auto exposure=diffLoan*rand()/RAND_MAX+minLoan;
		totalExposure+=exposure;
		return creditriskutilities::loan(getP(maxP), exposure, getNormalizedWeights(m));
	});
	lambda=totalExposure*lambda; //proxy for n*exposure*lambda
	q=q/lambda;

	const double xmax=0;
	const double xmin=-totalExposure*bL*maxP*.5*5;//5 is arbitrary

	auto expectation=creditriskutilities::computeExpectationVasicek(y0, alpha, tau);
	auto variance=creditriskutilities::computeVarianceVasicek(alpha, sigma, rho, tau);
	auto density=fangoost::computeInv(xNum, uNum,  xmin, xmax, [&](const auto& u){
		return creditriskutilities::executeVasicekMGF(
				creditriskutilities::logLPMCF(
					u,
					loans,
					[&](const auto& u, const auto& l){
						return creditriskutilities::lgdCF(u, l.exposure, alphL, bL, sigL, tau, bL);
					},
					lambda,
					q
				),
				expectation,
				variance
			);
	});
	auto dx=fangoost::computeDX(xNum, xmin, xmax);
	printJson(density, xmin, dx);

}
