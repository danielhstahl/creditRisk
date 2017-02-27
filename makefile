INCLUDES= -I ../FangOost -I ../FunctionalUtilities -I../rapidjson/include/rapidjson
creditRisk:main.o 
	g++ -std=c++14 -O3 -pthread main.o $(LDFLAGS) $(INCLUDES) -o creditRisk -fopenmp
main.o: main.cpp CreditUtilities.h
	g++ -std=c++14 -O3 -pthread  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o creditRisk
