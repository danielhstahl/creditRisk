
std::vector<Complex> logLPMCF(const Complex &u, const std::vector<loan>& loans, const auto& lgdCF, double lambda, double q){
	int n=loans.size();
	int m=loans[0].w.size();
	std::vector<Complex> v;//=std::vector<Complex>(m);
	Complex upperU=-lambda*u; //liquidity risk..u*lambda*i
	//upperU=upperU.exp().subtract(1.0).multiply(q);//.multiply(new Complex(0, 1));//q*(exp(i*u*lambda)-1)
	upperU=(exp(upperU)-1.0)*q;//.multiply(new Complex(0, 1));//q*(exp(i*u*lambda)-1)
	upperU=upperU-u; //liquidity..u*i+q*(exp(u*lambda*i)-1)
	std::vector<Complex> helperPhi;//=std::vector<Complex>(n);
	v.emplace_back(Complex(0, 0));
	for(int i=0; i<n;i++){
		helperPhi.emplace_back(lgdCF(upperU, loans[i])-1.0)*loans[i].pd);//e^{i*u*l_j}-1
		v[0]=v[0]+helperPhi[i]*loans[i].w[0]; //w_{j, k}p_j*(e^{i*u*l_j}-1)
	}
	for(int j=1; j<m; j++){
		v.emplace_back(Complex(0, 0));
		for(int i=0; i<n;i++){
			v[j]=v[j]+helperPhi[i]*loans[i].w[j]; //w_{j, k}p_j*(e^{i*u*l_j}-1)
		}
	}
	return(v);
}
