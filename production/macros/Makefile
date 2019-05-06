os = $(shell uname -s)

#INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(STARPICOPATH)
INCFLAGS      = -I$(shell root-config --incdir) -I$(FASTJETDIR)/include -I$(STARPICOPATH) -I$(ROOUNFOLDPATH) -I/opt/local/include

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

LIBPATH       = $(ROOTLIBS) $(FJLIBS) -L$(STARPICOPATH) -L$(ROOUNFOLDDIR) #-L$(FASTJETDIR)/lib -L$(STARPICOPATH) -L$(ROOUNFOLDDIR)
LIBS          = -lRecursiveTools -lfastjettools -lfastjet -lTStarJetPico -lRooUnfold


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
	$(CXX) $(LDFLAGS) $(LIBPATH) $^ $(LIBS) -o $@

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all : $(BDIR)/construct_response

response : $(BDIR)/construct_response

$(ODIR)/funcs.o         : $(SDIR)/funcs.cxx $(SDIR)/funcs.hh
$(ODIR)/construct_response.o		: $(SDIR)/construct_response.cxx

#data analysis
$(BDIR)/construct_response		: $(ODIR)/construct_response.o $(ODIR)/funcs.o
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

