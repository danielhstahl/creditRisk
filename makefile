LDFLAGS=-L ../Complex -lComplex -L ../FangOosterlee -lFangOosterlee
INCLUDES=-I ../Complex  -I ../FangOosterlee
creditRisk:main.o lgdCF.o 
	g++ -std=c++14 -O3  main.o lgdCF.o $(LDFLAGS) $(INCLUDES) -o creditRisk -fopenmp
main.o: main.cpp IntegroVasicekMG.h lgdCF.h lpmCF.h
	g++ -std=c++14 -O3  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
lgdCF.o: lgdCF.cpp lgdCF.h
	g++ -std=c++14 -O3 -c lgdCF.cpp $(LDFLAGS) $(INCLUDES) -fopenmp


clean:
	     -rm *.o creditRisk
