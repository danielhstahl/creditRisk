#include "lpmCF.h"

lpmCF::lpmCF(){
}

std::vector<Complex> lpmCF::logCF(Complex &u, std::vector<double > &p, std::vector<std::map<std::string, double> > &l, std::vector<std::vector<double> > &w, double lambda, double q){
	int n=p.size();
	int m=w[0].size();
	std::vector<Complex> v=std::vector<Complex>(m);
	Complex upperU=u.multiply(-lambda); //liquidity risk..u*lambda*i
	upperU=upperU.exp().subtract(1.0).multiply(q);//.multiply(new Complex(0, 1));//q*(exp(i*u*lambda)-1)
	upperU=upperU.subtract(u); //liquidity..u*i+q*(exp(u*lambda*i)-1)
	std::vector<Complex> helperPhi=std::vector<Complex>(n);
	v[0]=Complex(0, 0);
	for(int i=0; i<n;i++){
		//std::map<std::string, double> li=l[i];
		lgdCF phiL(l[i]["exposure"], l[i]["alpha"], l[i]["beta"], l[i]["sigma"], l[i]["tau"],l[i]["x0"]);
		//helperPhi[i]=Complex(w[i][0]*p[i], 0); //w_{j, k}*p_j
		helperPhi[i]=phiL.execute(upperU).subtract(1.0).multiply(p[i]);//e^{i*u*l_j}-1
		v[0]=v[0].add(helperPhi[i].multiply(w[i][0])); //w_{j, k}p_j*(e^{i*u*l_j}-1)
	}
	for(int j=1; j<m; j++){
		v[j]=Complex(0, 0);
		for(int i=0; i<n;i++){
			v[j]=v[j].add(helperPhi[i].multiply(w[i][j])); //w_{j, k}p_j*(e^{i*u*l_j}-1)
		}
	}

	//std::vector<Complex>().swap(helperPhi);//clear memory...but doesn't work?
	return(v);
}
/*int lpmCF::getM(){
	return m;
}*/
