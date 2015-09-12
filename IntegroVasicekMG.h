#ifndef __INTEGROVASICEKMG_H_INCLUDED__
#define __INTEGROVASICEKMG_H_INCLUDED__

#include <cmath>
#include <vector>
#include "Complex.h"
#include <iostream>


class IntegroVasicekMG {
	private:
		std::vector<double> alpha;
		std::vector<double> y0;
		std::vector<double > sigma;
		std::vector<std::vector<double> > rho;
		std::vector<std::vector<double> > variance;
		std::vector<double > expectation;
		double tau;
		int m;
		double helpComputeMoments(double);
		void computeMGF();
	public:
		IntegroVasicekMG();
		IntegroVasicekMG(std::vector<double>&, std::vector<double>&, std::vector<std::vector<double> >&, std::vector<double>& , double);
		Complex execute(const std::vector<Complex>&);
};
#endif
