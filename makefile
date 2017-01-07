INCLUDES= -I ../FangOost -I ../FunctionalUtilities -I../rapidjson/include/rapidjson
creditRisk:main.o 
	g++ -std=c++14 -O3  main.o $(LDFLAGS) $(INCLUDES) -o creditRisk -fopenmp
main.o: main.cpp IntegroVasicekMG.h lpmCF.h 
	g++ -std=c++14 -O3  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o creditRisk
