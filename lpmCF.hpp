
std::vector<Complex> logLPMCF(Complex &u, std::vector<loan>& loans, const auto& lgdCF, double lambda, double q){
	int n=loans.size();
	int m=loans[0].w.size();
	std::vector<Complex> v=std::vector<Complex>(m);
	Complex upperU=u.multiply(-lambda); //liquidity risk..u*lambda*i
	upperU=upperU.exp().subtract(1.0).multiply(q);//.multiply(new Complex(0, 1));//q*(exp(i*u*lambda)-1)
	upperU=upperU.subtract(u); //liquidity..u*i+q*(exp(u*lambda*i)-1)
	std::vector<Complex> helperPhi=std::vector<Complex>(n);
	v[0]=Complex(0, 0);
	for(int i=0; i<n;i++){
		helperPhi[i]=lgdCF(upperU, loans[i]).subtract(1.0).multiply(loans[i].pd);//e^{i*u*l_j}-1
		v[0]=v[0].add(helperPhi[i].multiply(loans[i].w[0])); //w_{j, k}p_j*(e^{i*u*l_j}-1)
	}
	for(int j=1; j<m; j++){
		v[j]=Complex(0, 0);
		for(int i=0; i<n;i++){
			v[j]=v[j].add(helperPhi[i].multiply(loans[i].w[j])); //w_{j, k}p_j*(e^{i*u*l_j}-1)
		}
	}
	return(v);
}
