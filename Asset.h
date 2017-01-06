#ifndef __ASSET_H_INCLUDED__
#define __ASSET_H_INCLUDED__

#include <cmath>
#include <vector>
#include "Complex.h"
struct loan{
    double pd;
    //std::function<Complex(const Complex&)> lgdCF;//characteristic function
    double exposure;
    std::vector<double> w;
    loan(double pd_, double exposure_, std::vector<double>&& w_){
        pd=pd_;
        exposure=exposure_;
        w=w_;
    }
};


#endif