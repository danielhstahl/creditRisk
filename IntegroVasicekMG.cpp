#include "IntegroVasicekMG.h"

IntegroVasicekMG::IntegroVasicekMG(){ //default constructor to ensure there is defaults for everything..
	m=1;
	tau=1;

}
IntegroVasicekMG::IntegroVasicekMG(std::vector<double> &alpha_, std::vector<double> &sigma_, std::vector<std::vector<double> > &rho_,std::vector<double> &y0_, double tau_) {
	alpha=alpha_;
	sigma=sigma_;
	rho=rho_;
	m=alpha.size();
	tau=tau_;
	y0=y0_;
}
double IntegroVasicekMG::helpComputeMoments(double alpha){ //hleper function since called so much
	return((1-std::exp(-alpha*tau))/alpha);
}

void IntegroVasicekMG::computeMGF(){

	variance.resize(m, std::vector<double>(m));
	for(int i=0; i<m; i++){
		double ai=helpComputeMoments(alpha[i]);
		expectation.push_back((y0[i]-1)*ai+tau);
		std::vector<double> row=std::vector<double>(m);
		for(int j=0; j<m; j++){
			double aj=helpComputeMoments(alpha[j]);
			double helpVarij=(rho[i][j]*sigma[i]*sigma[j]/(alpha[i]*alpha[j]))*(tau-ai-aj+(1-std::exp(-(alpha[i]+alpha[j])*tau))/(alpha[i]+alpha[j])); //difra page 10
			row[j]=helpVarij;
		}
		variance[i]=row;
	}

}

Complex IntegroVasicekMG::execute(const std::vector<Complex> &v){
	//int mv=v.length;
	Complex el=Complex(0, 0);
	Complex var=Complex(0, 0);
	if(variance.empty()){
		computeMGF();
	}
	for(int i=0; i<m; i++){
		el=el.add(v[i].multiply(expectation[i]));
		for(int j=0; j<m; j++){
			var=var.add(v[i].multiply(v[j]).multiply(variance[i][j]));
		}
	}
	var=var.multiply(.5).add(el);
	Complex phi=var.exp();
	//std::cout<<v[0].getReal()<<std::endl;
	return phi;
}
