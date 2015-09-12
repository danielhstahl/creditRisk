#ifndef __LPMCF_H_INCLUDED__
#define __LPMCF_H_INCLUDED__

//class lgdCF;
#include <cmath>
#include "Complex.h"
#include "lgdCF.h"
#include <vector>
#include <map>
#include <iostream> //debugging

class lpmCF { //implements ICharacteristicFunction
	private:
		/*std::vector<std::shared_ptr<double> > p; //"pd"
		std::vector<std::shared_ptr<std::map<std::string, double> > > l; //l is a parameter characterizing the LGD.  In most basic form, it is the dollar exposure
		std::vector<std::vector<std::shared_ptr<double> > > w; //weights for each asset with respect to the systemic variables
		std::vector<std::shared_ptr<double> > b; //balance on assets
		std::vector<std::shared_ptr<double> > r;//loss given liquidity
		double lambda0; //loss given liqudity inherent in the portfolio
		double lambda;
		double q; //"probability" of liquidity crisis
		//double tau;
		int n;
		int m;
		//double x0_;
		void init(std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<std::map<std::string, double> > >, std::vector<std::vector<std::shared_ptr<double> > > , std::vector<std::shared_ptr<double> > , std::vector<std::shared_ptr<double> > , double, double);*/
	public:
		lpmCF();
		/*lpmCF(std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<std::map<std::string, double> > >); //p, l
		lpmCF(std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<std::map<std::string, double> > >, double, double); //p, l, lambda0, q
		lpmCF(std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<std::map<std::string, double> > >, std::vector<std::vector<std::shared_ptr<double> > >, double, double); //p, l, w, lambda0, q
		lpmCF(std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<std::map<std::string, double> > >, std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<double> > , double, double); //p, l, b, r, lambda0, q,
		lpmCF(std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<std::map<std::string, double> > >, std::vector<std::vector<std::shared_ptr<double> > >); //p, l, w
		lpmCF(std::vector<std::shared_ptr<double> >, std::vector<std::shared_ptr<std::map<std::string, double> > >, std::vector<std::vector<std::shared_ptr<double> > > ,std::vector<std::shared_ptr<double> > , std::vector<std::shared_ptr<double> > , double, double); //p, l, w, b, r, lambda0, q
		std::vector<Complex> logCF(Complex);*/
		std::vector<Complex> logCF(Complex&, std::vector<double >&, std::vector<std::map<std::string, double> >&, std::vector<std::vector<double> >&, double, double);
		
		//int getM();
};
#endif