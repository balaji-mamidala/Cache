CXXFLAGS := -g -Wall -std=c++0x -lm
CXX=gcc

all: cachesim

cachesim: cachesim.o cachesim_driver.o
	$(CXX) -o cachesim cachesim.o cachesim_driver.o -lm

clean:
	rm -f cachesim *.o
