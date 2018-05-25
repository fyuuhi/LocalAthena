EXENAME = RPC

OBJS = RPC.o

AR = g++
CXX = g++

ROOTINCS      = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --glibs)

CXXFLAGS      = -g -Wall -fPIC $(ROOTINCS)
SOFLAGS       = -shared

all: $(OBJS)
	$(CXX) -o $(EXENAME) $(OBJS) $(ROOTLIBS)

%.o:%.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o:%.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o:%.cc
