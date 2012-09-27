.DELETE_ON_ERROR:
.PHONY: all clean checkdirs

ROOTCONFIG  := root-config
ROOTGLIBS   := $(shell $(ROOTCONFIG) --glibs) #-lcudart -lcuda
INCDIR      := $(shell $(ROOTCONFIG) --incdir)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs) -lSpectrum
CXXFLAGS    := -O2 $(ROOTCFLAGS) -m64#-arch=sm_12 -m64
SRC_OBJ     := $(addprefix $(OBJDIR)/,$(SRC_DIR:../.cu=.o))

INCLUDES    := -I $(INCDIR) 
NVCC        := g++

all: checkdirs plot

plot: plot.cc
	$(NVCC) $(CXXFLAGS) $(INCLUDES) $^ $1 -o $@ $(ROOTGLIBS) $(ROOTLIBS)
		        
clean:
	@rm -f *.o
	@rm -f fitNofit
