#include "lgdCF.h"

Complex lgdCF(const Complex &u, double &l, double &lambda, double &theta, double &sigma, double &t, double &x0){/*I think this is a CIR characteristic function*/
	double expt=exp(-lambda*t);
	double sigL=-sigma*sigma/(2*lambda);
	Complex uu=u*l;
	Complex uP=uu*(1-expt)*sigL+1.0;
	return exp((uu*expt*x0)/uP)*pow(uP, theta/sigL);
}
