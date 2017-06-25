INCLUDES= -I ../FangOost -I ../FunctionalUtilities  -I ../Vasicek -I ../rapidjson/include/rapidjson
creditRisk:main.o 
	g++ -std=c++14 -O3 $(STATIC) -pthread main.o $(INCLUDES) -o creditRisk -fopenmp
main.o: main.cpp CreditUtilities.h
	g++ -std=c++14 -O3 $(STATIC) -pthread  -c main.cpp  $(INCLUDES) -fopenmp
clean:
	     -rm *.o creditRisk
test: test.cpp CreditUtilities.h
	g++ -std=c++14 -pthread test.cpp $(INCLUDES) -o test -fopenmp