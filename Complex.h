#ifndef __COMPLEX_H_INCLUDED__
#define __COMPLEX_H_INCLUDED__

#include <cmath>
#include <iostream> //for debugging

class Complex {
	private:
		double real;
		double im;
	public:
		Complex();
		Complex(double, double);
		Complex multiply(double, double) const;
		Complex multiply(Complex) const;
		Complex multiply(double) const;
		Complex exp();
		Complex log();
		Complex add(Complex) const;
		Complex add(double) const;
		Complex subtract(Complex) const;
		Complex subtract(double) const;
		Complex power(double) const;
		double getReal();
		double getIm();
		Complex divide(Complex) const;
		Complex divide(double) const;
};

#endif
