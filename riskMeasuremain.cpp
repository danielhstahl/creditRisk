#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOosterleeNew.h"
#include "Complex.h"
#include <cmath>
#include <fstream>
#include "lpmCF.h"
#include "lgdCF.h"
#include "Asset.h"
#include <unordered_map>
#include "IntegroVasicekMG.h"
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono


//Complex LPMToInvert(Complex u,
int main(){

    int n=100000;
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
    double VaRAlpha=.01;//lower part
    computeVasicekMGF(vasicekEL, vasicekVR, alpha, sigma, rho, y0, tau);
    auto vkEL=[&](const double& u, const double& c, const double& d, const double& a){//d=d-a, so have to adjust
      if(u==0){
        return ((d+a)*(d+a)-c*c)*.5;
      }
      else{
        double intern1=u*d;
        double intern2=u*(c-a);
        return ((d+a)*sin(intern1)-c*sin(intern2))/u+(cos(intern1)-cos(intern2))/(u*u);
      }
    };
    auto vkMinusEL=[&](const double& u, const double& c, const double& d, const double& a){//d=d-a, so have to adjust
      if(u==0){
        return -((d+a)*(d+a)-c*c)*.5;
      }
      else{
        double intern1=u*d;
        double intern2=u*(c-a);
        return -((d+a)*sin(intern1)-c*sin(intern2))/u-(cos(intern1)-cos(intern2))/(u*u);
      }
    };
    auto vkCDF=[&](const double& u, const double& c, const double& d, const double& a){ //u=kpi/(b-a)
      if(u==0){
        return (d+a-c);//was d-c but already normalized
      }
      else{
        double intern1=u*d;//was d-a but already normalized
        double intern2=u*(c-a);
        return (sin(intern1)-sin(intern2))/u;
      }
    };
    std::ofstream outputCSV;
    outputCSV.open("creditDistribution.csv");
    outputCSV <<"x, y"<<std::endl;
    auto callbackEL=[&](const std::vector<double>& f, const double& du, const double& dx, const int& k, const int& h, const double& a){
      std::vector<double> holdConv(h);
      //double myExpectation=0;
      for(int i=0; i<h; ++i){
        holdConv[i]=0;
        for(int j=0; j<k; ++j){
          holdConv[i]+=vkEL(du*j, a, dx*i, a)*f[j];//*cos(du*j*dx*i);
        }
        outputCSV<<a+dx*i<<", "<<holdConv[i]<<std::endl;
      }
      outputCSV.close();
    };
    auto callbackCDF=[&](std::vector<double>& f, const double& du, const double& dx, const int& k, const int& h, const double& a){
      std::vector<double> holdConv(h);
      for(int i=0; i<h; ++i){
        holdConv[i]=0;
        for(int j=0; j<k; ++j){
          holdConv[i]+=vkCDF(du*j, a, dx*i, a)*f[j];//*cos(du*j*dx*i);//*cos(du*j*a);
          //holdConv[i]+=f[j]*(cos(du*j*dx*i));//*cos(du*j*a);;
        }
        outputCSV<<a+dx*i<<", "<<holdConv[i]<<std::endl;
        //std::cout<<"cdf "<<holdConv[i]<<std::endl;
      }
      outputCSV.close();
    };
    auto callbackExpectedShortfall=[&](std::vector<double>& f, const double& du, const double& dx, const int& k, const int& h, const double& a){
      double cdf=0;
      int i=0;
      while(cdf<VaRAlpha){
        cdf=0;
        for(int j=0; j<k; ++j){
          cdf+=vkCDF(du*j, a, dx*i, a)*f[j];

          //holdConv[i]+=f[j]*(cos(du*j*dx*i));//*cos(du*j*a);;
        }
        i++;
        //outputCSV<<a+dx*i<<", "<<holdConv[i]<<std::endl;
        //std::cout<<"cdf "<<holdConv[i]<<std::endl;
      }
      double expectedShortfall=0;
      for(int j=0; j<k; ++j){
        expectedShortfall+=vkMinusEL(du*j, a, dx*(i-1), a)*f[j];

        //holdConv[i]+=f[j]*(cos(du*j*dx*i));//*cos(du*j*a);;
      }
      std::cout<<expectedShortfall/VaRAlpha<<std::endl;
      //outputCSV.close();
    };
    FangOost invert(xmin, xmax, 128, 1024);
    invert.computeInv(
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
    //invert.computeExpectation(callbackCDF);
    invert.computeExpectation(callbackExpectedShortfall);
	auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	/*std::cout<<"actualEL: "<<invert.getEL()<<std::endl;
	std::ofstream outputCSV;
	outputCSV.open("creditDistribution2.csv");
	outputCSV <<"x, y"<<std::endl;
	for(int i=0; i<1024; i++){
		outputCSV<<""<<results["x"][i]<<", "<<results["y"][i]<<std::endl;
	}
	outputCSV.close();*/

}
