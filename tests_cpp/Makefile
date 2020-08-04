CC=g++
FLAGS=-g -std=c++17 -I/home/jhk/anaconda3/include/python3.7m/ -fopenmp
LIBS=-lfftw3 -lm -lgomp
GTEST_FLAGS=`pkg-config --cflags gtest_main`
GTEST_LIBS=`pkg-config --libs gtest_main`
CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(patsubst %.cpp, %.o, $(CPP_FILES))

.PHONY: all tests clean

all: testapp

tests: all
	./testapp

testapp: $(OBJ_FILES) ../src/Integrators.o ../src/kde.o
	$(CXX) $(FlAGS) $(GTEST_FLAGS) $^ -o $@ $(LIBS) $(GTEST_LIBS)

%.o: %.cpp
	$(CXX) $(FLAGS) $(GTEST_FLAGS) $< -c -o $@


clean:
	-rm ../src/*.o *.o testapp