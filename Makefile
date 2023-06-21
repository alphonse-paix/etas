CXX = g++
CXXFLAGS = -std=c++17 -Ofast -Wall -Wextra

etas: main.cpp
	$(CXX) $(CXXFLAGS) -I/home/alphonse/Documents/include main.cpp -o etas

