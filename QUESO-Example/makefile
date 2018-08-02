#module restore cheese
#export LD_LIBRARY_PATH=/opt/apps/ossw/libraries/queso/queso-0.55.0/c7/gcc-5.2/mpich2-3.1.4/boost-1.58.0/gsl-2.0/hdf5-1.8.15-patch1/lib:\
#/opt/apps/ossw/libraries/boost/boost-1.58.0/c7/gcc-5.2/lib:\
#/opt/apps/ossw/libraries/gsl/gsl-2.0/c7/gcc-5.2/lib:\
#/opt/apps/ossw/libraries/hdf5/hdf5-1.8.15-patch1/c7/gcc-5.2/lib
#make clean
#make -f makefile
QUESO_DIR = /opt/apps/ossw/libraries/queso/queso-0.55.0/c7/gcc-5.2/mpich2-3.1.4/boost-1.58.0/gsl-2.0/hdf5-1.8.15-patch1
BOOST_DIR = /opt/apps/ossw/libraries/boost/boost-1.58.0/c7/gcc-5.2
GSL_DIR = /opt/apps/ossw/libraries/gsl/gsl-2.0/c7/gcc-5.2
HDF5_DIR = /opt/apps/ossw/libraries/hdf5/hdf5-1.8.15-patch1/c7/gcc-5.2

INC_PATHS = \
	-I. \
	-I$(QUESO_DIR)/include \
	-I$(BOOST_DIR)/include \
	-I$(GSL_DIR)/include \
	-I$(HDF5_DIR)/include

LIBS = \
	-L$(QUESO_DIR)/lib -lqueso \
	-L$(BOOST_DIR)/lib -lboost_program_options \
	-L$(GSL_DIR)/lib -lgsl \
	-L$(HDF5_DIR)/lib -lhdf5

CXX = mpic++
CXXFLAGS += -g -Wall -c #-std=c++11

default: all

.SUFFIXES: .o .C

all:	gsl

clean:
	rm -f *~
	rm -f *.o
	rm -f gsl
	rm -f outputData/*
	rm -f *.jpg

gsl: main.o likelihood.o compute.o
	$(CXX) main.o \
	likelihood.o \
	compute.o \
	-o gsl $(LIBS)

%.o: %.C
	$(CXX) $(INC_PATHS) $(CXXFLAGS) $<
