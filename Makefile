os = $(shell uname -s)

#INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(STARPICOPATH)
INCFLAGS      = -I$(shell root-config --incdir) -I$(FASTJETDIR)/include -I$(STARPICOPATH) -I/opt/local/include

ifeq ($(os),Linux)
CXXFLAGS      = -std=c++11
else
CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init -std=c++11
## for debugging:
# CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared 
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = clang
endif


ROOTLIBS      = $(shell root-config --libs)
FJLIBS	      = $(shell fastjet-config --libs)
#PYTHIALIBS    = $(shell pythia8-config --ldflags)

LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(STARPICOPATH)
LIBS          = -lfastjet -lfastjettools -lTStarJetPico


# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################


###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all : $(BDIR)/ppsim $(BDIR)/ppdata 

data : $(BDIR)/ppdata
sim : $(BDIR)/ppsim

#$(SDIR)/dict.cxx                : $(SDIR)/ktTrackEff.hh
#	cd ${SDIR}; rootcint -f dict.cxx -c -I. ./ktTrackEff.hh

#$(ODIR)/dict.o                  : $(SDIR)/dict.cxx
#$(ODIR)/ktTrackEff.o            : $(SDIR)/ktTrackEff.cxx $(SDIR)/ktTrackEff.hh

$(ODIR)/funcs.o		: $(SDIR)/funcs.cxx $(SDIR)/funcs.hh
$(ODIR)/ppsim.o		: $(SDIR)/ppsim.cxx
$(ODIR)/ppdata.o	: $(SDIR)/ppdata.cxx

#data analysis
$(BDIR)/ppsim		: $(ODIR)/ppsim.o $(ODIR)/funcs.o #$(ODIR)/ktTrackEff.o $(ODIR)/dict.o
$(BDIR)/ppdata		: $(ODIR)/ppdata.o $(ODIR)/funcs.o

###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*
	rm -fr log/*

