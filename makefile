LDFLAGS=-L ../Complex -lComplex -L ../FangOosterlee -lFangOosterlee
INCLUDES=-I ../Complex  -I ../FangOosterlee
creditRisk:main.o lgdCF.o lpmCF.o IntegroVasicekMG.o
	g++ -std=c++11 -O3  main.o lgdCF.o lpmCF.o IntegroVasicekMG.o  $(LDFLAGS) $(INCLUDES) -o creditRisk -fopenmp
main.o: main.cpp IntegroVasicekMG.h lgdCF.h lpmCF.h
	g++ -std=c++11 -O3  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
lgdCF.o: lgdCF.cpp lgdCF.h
	g++ -std=c++11 -O3 -c lgdCF.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
lpmCF.o: lpmCF.cpp lpmCF.h
	g++ -std=c++11 -O3 -c lpmCF.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
IntegroVasicekMG.o: IntegroVasicekMG.cpp IntegroVasicekMG.h
	g++ -std=c++11 -O3 -c IntegroVasicekMG.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o creditRisk
