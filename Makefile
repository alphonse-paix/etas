CXX = g++
CXXFLAGS = -std=c++17 -Ofast -Wall -Wextra
INCLPATH = ~/Documents/include

etas: main.cpp
	$(CXX) $(CXXFLAGS) -I $(INCLPATH) main.cpp -o etas
