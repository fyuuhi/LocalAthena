EXENAME = RPC

OBJS = RPC.o

AR = g++
CXX = g++

ROOTINCS      = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --glibs)

CXXFLAGS      = -g -Wall -fPIC $(ROOTINCS)
SOFLAGS       = -shared

all: $(OBJS)
	$(CXX) -I /gpfs/home/yfukuhar/tmp/plog/include -o $(EXENAME) $(OBJS) $(ROOTLIBS)

%.o:%.cxx
	$(CXX) $(CXXFLAGS) -I /gpfs/home/yfukuhar/tmp/plog/include -c $< -o $@

%.o:%.C
	$(CXX) $(CXXFLAGS) -I /gpfs/home/yfukuhar/tmp/plog/include -c $< -o $@

%.o:%.cc

clean:
	rm -f $(EXENAME) $(OBJS)
