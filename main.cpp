#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOosterlee.h"
#include "Complex.h"
#include <cmath>
#include "lpmCF.h"
#include "lgdCF.h"
#include "Asset.h"
#include <unordered_map>
#include "IntegroVasicekMG.h"
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson

int main(){
	int xNum=1024;
	int uNum=256; 
  int n=100000;//sensible defaults
	int m=1;
	std::vector<double> alpha=std::vector<double>(m);
	std::vector<double> sigma=std::vector<double>(m);
	std::vector<double> y0=std::vector<double>(m);
	std::vector< std::vector<double> > rho(m, std::vector<double> (m));
	double tau=1;
	alpha[0]=.1;
	sigma[0]=.3;
	rho[0][0]=1;
	y0[0]=1;

	/*	sigma[0]=.2;
		sigma[1]=.1;
		sigma[2]=.3;*/
	/*alpha[0]=.2;
	alpha[1]=.3;
	alpha[2]=.4;

	sigma[0]=.6;
	sigma[1]=.6;
	sigma[2]=.6;
	y0[0]=.9;
	y0[1]=1;
	y0[2]=1.1;
	rho[0][0]=1;
	rho[1][1]=1;
	rho[2][2]=1;
	rho[0][1]=.2;
	rho[1][0]=.2;
	rho[0][2]=-.3;
	rho[2][0]=-.3;
	rho[1][2]=.1;
	rho[2][1]=.1;*/


	//These are parameters for the LGD CF function
  double alphL=.2;
	double bL=.5;
	double sigL=.2;
	double sumL=0;
	srand(5);
	double maxP=.09;
	double el=0;
  double totalwP=0;
	double q=.05;
	double lambda=.05;
	auto runParameters=[&](std::string& parameters){
		rapidjson::Document parms;
		parms.Parse(parameters.c_str());//yield data
		parameters.clear();
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
		std::vector<loan> loans;
		for(int i=0; i<n; i++){
	    double exposure=40000.0*rand()/RAND_MAX+10000.0;
	    std::vector<double> w(m);
	    totalwP=0;
	    for(int j=0; j<m; j++){
				w[j]=1.0*rand(); //randomly assign weights
				totalwP=totalwP+w[j];
			}
			for(int j=0; j<m; j++){
				w[j]=w[j]/totalwP;
			}
	    sumL+=exposure;
	    loans.push_back(loan(maxP*rand()/RAND_MAX+.0001,exposure,w));//probability constrained between .0001 and maxP
			el=el+exposure*loans[i].pd*bL;
		}
		FangOosterlee invert(uNum, xNum);
		lambda=lambda*sumL;
		q=q/sumL;
		double xmax=0;
		double xmin=-sumL*bL*maxP*.5*5;//5 is arbitrary
		/*std::cout<<"EL: "<<-el<<std::endl;*/
	  std::vector<double> vasicekEL;
	  std::vector<std::vector<double> > vasicekVR;
	  computeVasicekMGF(vasicekEL, vasicekVR, alpha, sigma, rho, y0, tau);
		invert.computeDistributionJSON(
			true,
      xmin,
      xmax,
	    [&](Complex& u) {
	      return executeVasicekMGF(
	        logLPMCF(
	          u,
	          loans,
	          [&](Complex& u, loan& l){
	            return lgdCF(u, l.exposure, alphL, bL, sigL, tau, bL);
	          },
	          lambda,
	          q
	        ),
	        vasicekEL,
	        vasicekVR
	      );
	    }
	  );
	};
	while(true){
		std::string parameters;
		for (parameters; std::getline(std::cin, parameters);) {
			break;
		}
		runParameters(parameters);
	}


}
