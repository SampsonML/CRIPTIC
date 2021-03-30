# Get location of AMREX library
AMREXLIBDIR := $(AMREX_LIBRARY_HOME)/lib
AMREXINCDIR := $(AMREX_LIBRARY_HOME)/include

# Get location of GSL library
GSLLIBDIR   := $(GSL_LIBRARY_HOME)/lib
GSLINCDIR   := $(GSL_LIBRARY_HOME)/include

# Read the flags that were used to compile AMRex
COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(AMREXLIBDIR)/pkgconfig/amrex.pc)
COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(AMREXLIBDIR)/pkgconfig/amrex.pc)

# Get location of HDF5 library
HDF5LIBDIR     := /home/krumholz/hdf5-1.8.21/hdf5/lib
HDF5INCDIR     := /home/krumholz/hdf5-1.8.21/hdf5/include

# Set compiler and compiler flags
CXX = mpicxx

CXXFLAGS := -I$(AMREXINCDIR) -I$(GSLINCDIR) -I$(HDF5INCDIR) $(COMPILE_CPP_FLAGS) 
LFLAGS := -lgomp -lgsl -lgslcblas -lhdf5 -lhdf5_cpp  -L$(AMREXLIBDIR) -L$(GSLLIBDIR) -L$(HDF5LIBDIR) $(COMPILE_LIB_FLAGS) 

# Pointers to sources and objects
SOURCES         = $(filter-out main.cpp, $(wildcard *.cpp))
OBJECTS         = $(SOURCES:%.cpp=%.o)
DEPS            = $(SOURCES:%.cpp=%.d)

# Default target
.PHONY: all criptic clean

all: criptic

# Include dependencies
-include $(DEPS)

criptic: $(OBJECTS) main.o
	$(CXX) -o criptic $^ $(LFLAGS)

clean:
	rm -f criptic $(OBJECTS) main.o Prob.cpp Definitions.H
