#include "lgdCF.h"
lgdCF::lgdCF(){

}
/*lgdCF::lgdCF(double l_, double alpha, double b, double sigma_, double t_, double x0_) { //dX=alpha(beta-x)dt+sigma sqrt(x) dW
	l=l_;
	lambda=alpha; //
	theta=b;
	//lambda=lambda_;
	sigma=sigma_;
	t=t_;
	x0=x0_;
}*/
Complex lgdCF::execute(Complex &u, double &l, double &lambda, double &theta, double &sigma, double &t, double &x0){
	double expt=std::exp(-lambda*t);
	double sigL=-sigma*sigma/(2*lambda);
	Complex uu=u.multiply(l);
	Complex uP=uu.multiply(1-expt).multiply(sigL).add(1.0);
	return uu.multiply(expt*x0).divide(uP).exp().multiply(uP.power(theta/sigL));


}
