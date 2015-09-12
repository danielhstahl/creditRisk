#include "Complex.h"

Complex::Complex() {
	real=0;
	im=0;
}
Complex::Complex(double x, double y) {
	real=x;
	im=y;
}
Complex Complex::multiply(double x, double y) const { //consider removing
	Complex plac=Complex(real*x-im*y, im*x+real*y);
	return plac;
}
Complex Complex::multiply(Complex c) const {
	Complex plac=Complex(real*c.getReal()-im*c.getIm(), im*c.getReal()+real*c.getIm());
	return plac;
}
Complex Complex::multiply(double c) const {
	Complex plac=Complex(c*real, c*im);
	return plac;
}
Complex Complex::exp() {
	Complex plac=Complex(std::exp(real)*std::cos(im), std::exp(real)*std::sin(im));
	return plac;
}
Complex Complex::log() {
	double modulus=std::sqrt(real*real+im*im);
	Complex plac=Complex(std::log(modulus), std::atan2(im, real));
	return plac;
}
Complex Complex::add(Complex c) const {
	Complex plac=Complex(real+c.getReal(), im+c.getIm());
	return plac;
}
Complex Complex::add(double c) const {
	Complex plac=Complex(real+c, im);
	return plac;
}
Complex Complex::subtract(double c) const {
	Complex plac=Complex(real-c, im);
	return plac;
}
Complex Complex::subtract(Complex c) const {
	Complex plac=Complex(real-c.getReal(), im-c.getIm());
	return plac;
}
Complex Complex::divide(Complex c) const {
	double den=c.getReal()*c.getReal()+c.getIm()*c.getIm();
	Complex plac=Complex((real*c.getReal()+im*c.getIm())/den, (im*c.getReal()-real*c.getIm())/den);
	return plac;
}
Complex Complex::divide(double c) const {
	Complex plac=Complex(real/c, im/c);
	return plac;
}
Complex Complex::power(double exponent) const {
	double modulus=std::sqrt(real*real+im*im);
	double arg=std::atan2(im, real);
	double log_re=std::log(modulus);
	double log_im=arg;
	double x_log_re=exponent*log_re;
	double x_log_im=exponent*log_im;
	double modulus_ans=std::exp(x_log_re);
	Complex plac=Complex(modulus_ans*std::cos(x_log_im), modulus_ans*std::sin(x_log_im));
	return plac;

}
double Complex::getReal() {
	return real;
}
double Complex::getIm() {
	return im;
}
