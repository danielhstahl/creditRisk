#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOosterlee.h"
#include "Complex.h"
#include <cmath>
#include <fstream>
#include "lpmCF.h"
#include "lgdCF.h"
#include "Asset.h"
#include "IntegroVasicekMG.h"
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono


//Complex LPMToInvert(Complex u,
int main(){
	FangOosterlee invert(256, 1024);
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
/*	sigma[0]=.2;
	sigma[1]=.1;
	sigma[2]=.3;*/
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
	rho[2][1]=.1;


	//These are parameters for the LGD CF function
    double alphL=.2;
	double bL=.5;
	double sigL=.2;
    
    std::vector<loan> loans;
	double sumL=0;
	srand(5);
	double maxP=.09;
	double el=0;
    double totalwP=0;
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
        loans.push_back(
            loan(
                maxP*rand()/RAND_MAX+.0001,//probability constrained between .0001 and maxP
                exposure,
                w
            )
        );
		el=el+exposure*loans[i].pd*bL;
	}
	double lambda=.05*sumL;
	double q=2.0/sumL;
	//lpmCF lp;
    //IntegroVasicekMG mg(alpha, sigma, rho, y0, tau);

	double xmax=0;
	double xmin=-sumL*bL*maxP*.5*5;//5 is arbitrary
	auto start = std::chrono::system_clock::now();
	std::cout<<"EL: "<<-el<<std::endl;
    
    //std::cout<<"Test for lgdCF: "<<loans[0].lgdCF(Complex(.5, .5)).getReal()<<std::endl;
    //double exposure=20989.8;
    //std::cout<<"Test for lgdCF: "<<lgdCF(Complex(.5, .5), exposure, alphL, bL, sigL, tau, bL).getReal()<<std::endl;
      
    std::vector<double> vasicekEL;
    std::vector<std::vector<double> > vasicekVR;
    computeVasicekMGF(vasicekEL, vasicekVR, alpha, sigma, rho, y0, tau);
	std::unordered_map<std::string, std::vector<double> > results=invert.computeDistribution(
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
	auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	std::cout<<"actualEL: "<<invert.getEL()<<std::endl;
	std::ofstream outputCSV;
	outputCSV.open("creditDistribution2.csv");
	outputCSV <<"x, y"<<std::endl;
	for(int i=0; i<1024; i++){
		outputCSV<<""<<results["x"][i]<<", "<<results["y"][i]<<std::endl;
	}
	outputCSV.close();

}
