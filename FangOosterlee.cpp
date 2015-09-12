#include "FangOosterlee.h"

FangOosterlee::FangOosterlee(int k_, int h_) {
	k=k_; //u discretions
	h=h_; //x discretions
	//M_PI=3.14159265358979323846;
}
double FangOosterlee::getEL() {
	return exloss;
}
double FangOosterlee::getVariance(){
	return vloss;
}
