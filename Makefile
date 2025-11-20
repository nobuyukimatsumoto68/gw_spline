CXX=g++
CXXFLAGS=-Wall -Wextra -std=c++23 -march=native -O3
INCLUDES=-I"/mnt/hdd_barracuda/opt/eigen/" -I"/mnt/hdd_barracuda/opt/HighFive/include/" -I"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/include/" -I"~/opt/eigen/"
LDFLAGS=-L"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/lib/" -L"/usr/lib/" -L"/usr/local/lib/" -lhdf5 -lgsl -lgslcblas -lm # -lhdf5_cpp
# INCLUDES=-I"/Users/nobuyukimatsumoto/opt/eigen/"
# LDFLAGS=

SRCS := $(wildcard *.cc)
OBJS := $(SRCS:%.cc=%.o)

all: $(OBJS)

%.d: %.cc
	$(CXX) $(INCLUDES) -M $< -o $@

include $(SRCS:.cc=.d)

%.o: %.cc
	$(CXX) $< -o $@ $(INCLUDES) $(LDFLAGS)

# all:
# 	g++ hist_spline.cc ${INCLUDES} ${LDFLAGS}
# 	g++ get_potential.cc ${INCLUDES} ${LDFLAGS}
