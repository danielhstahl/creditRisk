auto helpComputeMoments(const auto& alpha, const auto& tau){ //hleper function since called so much
	return (1-exp(-alpha*tau))/alpha;
}
template<typename number>
void computeVasicekMGF(std::vector<number> &expectation, std::vector< std::vector<number> > &variance, const std::vector<double> &alpha, const std::vector<double> &sigma, const std::vector<std::vector<double> > &rho, const std::vector<double> &y0, double tau){
	int m=alpha.size();
    for(int i=0; i<m; i++){
		auto ai=helpComputeMoments(alpha[i], tau);
		expectation.emplace_back((y0[i]-1)*ai+tau);
		std::vector<double> row=std::vector<double>(m);
		for(int j=0; j<m; j++){
			double aj=helpComputeMoments(alpha[j], tau);
			double helpVarij=(rho[i][j]*sigma[i]*sigma[j]/(alpha[i]*alpha[j]))*(tau-ai-aj+(1-std::exp(-(alpha[i]+alpha[j])*tau))/(alpha[i]+alpha[j])); //difra page 10
			row[j]=helpVarij;
		}
		variance.emplace_back(row);
	}
}
template<typename number>
Complex executeVasicekMGF(const std::vector<Complex> &v, const std::vector<number>& expectation , const std::vector< std::vector<number> >& variance){
	Complex el=Complex(0, 0);
	Complex var=Complex(0, 0);
    int m=expectation.size();
	for(int i=0; i<m; i++){
		el=el+v[i]*expectation[i];
		for(int j=0; j<m; j++){
			//var=var.add(v[i].multiply(v[j]).multiply(variance[i][j]));
			var=var+v[i]*v[j]*variance[i][j];
		}
	}
	var=var*.5+el;
	//Complex phi=var.exp();
	return exp(var);
}
