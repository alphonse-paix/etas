CXX = g++
CXXFLAGS = -std=c++17 -Ofast -Wall -Wextra
INCLPATH = ~/Include

etas: main.cpp
	$(CXX) $(CXXFLAGS) -I $(INCLPATH) main.cpp -o etas
