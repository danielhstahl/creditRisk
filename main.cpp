#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOosterlee.h"
#include "Complex.h"
#include <cmath>
#include <fstream>
#include "lpmCF.h"
#include "lgdCF.h"
#include "IntegroVasicekMG.h"
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono


//Complex LPMToInvert(Complex u,
int main(){
	FangOosterlee invert(128, 1024);
  int n=10000000;
	int m=3;
	std::vector<double> alpha=std::vector<double>(m);
	std::vector<double> sigma=std::vector<double>(m);
	std::vector<double> y0=std::vector<double>(m);
	std::vector< std::vector<double> > rho(m, std::vector<double> (m));
	double tau=1;
	alpha[0]=.2;
	alpha[1]=.3;
	alpha[2]=.4;
	sigma[0]=.15;
	sigma[1]=.1;
	sigma[2]=.2;
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
	rho[2][1]=.1;


	//double x0=.6;
	//IntegroVasicekMG systemicRisk=IntegroVasicekMG(alpha, sigma, rho, y0, tau);

	std::vector<double> p=std::vector<double>(n);
	std::vector<std::map<std::string, double> > l=std::vector<std::map<std::string, double> >(n);
	double alphL=.2;
	double bL=.5;
	double sigL=.2;
	std::vector<std::vector<double > > w(n, std::vector<double >(m) );
	double sumL=0;
	srand(5);
	double maxP=.09;
	double el=0;
	for(int i=0; i<n; i++){
		p[i]=maxP*rand()/RAND_MAX+.0001; //probability constrained between .0001 and .0301
		//delete pd;
		std::map<std::string, double> lMap;
		lMap["exposure"]=40000.0*rand()/RAND_MAX+10000.0; //exposures between 10000 and 50000...remember that these are just exposures, there is additional variability around this via the CIR process.  This value scales the CIR process multiplicatively
		lMap["alpha"]=alphL;
		lMap["beta"]=bL;
		lMap["sigma"]=sigL;
		lMap["tau"]=tau;
		lMap["x0"]=bL;
		l[i]=lMap;
		sumL=lMap["exposure"]+sumL;
		std::vector<double > wP=std::vector<double>(m);
		double totalwP=0;
		std::vector<double> wpj(m, 0);
		for(int j=0; j<m; j++){
			wpj[j]=1.0*rand(); //randomly assign weights
			totalwP=totalwP+wpj[j];
		}
		for(int j=0; j<m; j++){
			double wj=wpj[j]/totalwP;
			wP[j]=wj;

		}
		w[i]=wP;
		el=el+lMap["exposure"]*p[i]*bL;
	}
	double lambda=.05*sumL;
	double q=2.0/sumL;
	//double lambda=0;
	//double q=0;

	lpmCF lp;
  IntegroVasicekMG mg(alpha, sigma, rho, y0, tau);

	double xmax=0;
	double xmin=-sumL*bL*maxP*.5*5;//5 is arbitrary
	auto start = std::chrono::system_clock::now();
	std::cout<<"EL: "<<-el<<std::endl;

	std::map<std::string, std::vector<double> > results=invert.computeDistribution(xmin, xmax, [&](Complex u) {return mg.execute(lp.logCF(u, p, l, w, lambda, q));});
	auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	std::cout<<"actualEL: "<<invert.getEL()<<std::endl;
	std::ofstream outputCSV;
	outputCSV.open("creditDistribution.csv");
	outputCSV <<"x, y"<<std::endl;
	for(int i=0; i<1024; i++){
		outputCSV<<""<<results["x"][i]<<", "<<results["y"][i]<<std::endl;
	}
	outputCSV.close();

}
