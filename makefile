#------------------ compiler --------------

## GNU C++ compiler
CXX:= g++
CXXFLAGS:=-std=c++11 -fopenmp -O3

## Intel compiler
# CXX              :=  icpc
# CXXFLAGS         :=  -O3 -qopenmp -Wfatal-errors

#paths of source files
DIR_PLY:=./src/Forward_polyhedron
DIR_TET:=./src/Forward_tetrahedron
DIR_GEOMETRY:=./src/Geometry
DIR:=./src

#include
INCLUDE_PLY:=-I$(DIR_PLY)
INCLUDE_TET:=-I$(DIR_TET)
INCLUDE_GEOMETRY:=-I$(DIR_GEOMETRY)
INCLUDE_OTHERS:=-I$(DIR)

#source files
SRCS:=$(wildcard $(DIR_PLY)/*.cpp) $(wildcard $(DIR_GEOMETRY)/*.cpp) $(wildcard $(DIR_TET)/*.cpp)
OBJS             := $(patsubst %.cpp, %.o, $(SRCS))
#---------------------------------------------------------------------------------------
all: GraPly GraTet Sites
.PHONY: all

#link object files to executable files
GraPly: $(OBJS) ./src/GraPly.o
	$(CXX) $(CXXFLAGS) ./src/GraPly.o $(OBJS) -o GraPly 
GraTet: $(OBJS) ./src/GraTet.o
	$(CXX) $(CXXFLAGS) ./src/GraTet.o $(OBJS) -o GraTet
Sites: ./src/sites_generation.cpp
	$(CXX) $(CXXFLAGS) ./src/sites_generation.cpp $(OBJS) -o Sites

#compile c++ files
%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCLUDE_GEOMETRY) $(INCLUDE_OTHERS) $(INCLUDE_PLY) $(INCLUDE_TET) -c $< -o $@
./src/GraPly.o: ./src/GraPly.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_GEOMETRY) $(INCLUDE_OTHERS) $(INCLUDE_PLY) ./src/GraPly.cpp -c -o ./src/GraPly.o
./src/GraTet.o: ./src/GraTet.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_GEOMETRY) $(INCLUDE_OTHERS) $(INCLUDE_PLY) $(INCLUDE_TET) ./src/GraTet.cpp -c -o ./src/GraTet.o

echo:	
	@echo -e "Source Files:\t$(SRCS)"
	@echo -e "Object Files:\t$(OBJS)"
	@echo -e "C++ Compiler:\t$(CXX)"
	@echo -e "CXXFLAGS:    \t$(CXXFLAGS)"	

.PHONY: clean
clean:
	@rm -rf *.o *~  $(OBJS) ./src/GraPly.o ./src/GraTet.o ./src/sites_generation.o GraPly GraTet Sites