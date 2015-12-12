#include "lgdCF.h"

Complex lgdCF(const Complex &u, double &l, double &lambda, double &theta, double &sigma, double &t, double &x0){/*I think this is a CIR characteristic function*/
	double expt=std::exp(-lambda*t);
	double sigL=-sigma*sigma/(2*lambda);
	Complex uu=u.multiply(l);
	Complex uP=uu.multiply(1-expt).multiply(sigL).add(1.0);
	return uu.multiply(expt*x0).divide(uP).exp().multiply(uP.pow(theta/sigL));
}
