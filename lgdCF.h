#ifndef __LGDCF_H_INCLUDED__
#define __LGDCF_H_INCLUDED__

#include <cmath>
#include <iostream> //for debugging
#include "Complex.h"

class lgdCF { 
	private:
		double l;
		double theta;
		double lambda;
		double sigma;
		double t;
		double x0;
	public:
		lgdCF(double, double, double, double, double, double);
		Complex execute(Complex&);
};
#endif