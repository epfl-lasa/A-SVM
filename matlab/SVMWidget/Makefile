# This Makefile is used under Linux


# for Mac
# MATLABDIR ?= /opt/local/matlab

CXX = g++
CFLAGS = -Wall  -O3 -fPIC -I$(MATLABDIR)/extern/include -I. -I../../include -fexceptions
LDFLAGS= $(CXXFLAGS)  -Wl,--rpath -Wl,$(shell pwd)/../../lib
MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)" -cxx -O LDOPTIMFLAGS="$(LDFLAGS)"
# comment the following line if you use MATLAB on 32-bit computer
MEX_OPTION += -largeArrayDims -L../../lib -lnlopt -lASVMLearning 
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)



all:	smo nlopt  

nlopt:	mxNLOPTSolver 

smo: mxSMOSolver 

mxSMOSolver:     mxSMOSolver.c  mxUtil.o
	$(MEX) $(MEX_OPTION) mxSMOSolver.c  mxUtil.o

mxNLOPTSolver:     mxNLOPTSolver.c  mxUtil.o
	$(MEX) $(MEX_OPTION) mxNLOPTSolver.c  mxUtil.o

mxUtil.o:  mxUtil.c  mxUtil.h
	$(CXX) $(CFLAGS) -c  mxUtil.c

clean:
	rm -f *~ *.o *.mex* *.obj
	
	