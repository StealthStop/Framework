WDIR = .
ODIR = obj
IDIR = $(CMSSW_BASE)/src/Framework/Framework/include
SDIR = $(CMSSW_BASE)/src/Framework/Framework/src
TDIR = $(CMSSW_BASE)/src

OPENCV_DIRECTORY = /uscms_data/d3/nstrobbe/SOFTWARE/opencv
TTDIR            = $(CMSSW_BASE)/src/TopTagger/
TTTDIR           = $(TTDIR)/TopTagger/test

CXX = g++

CXXFLAGS  += -I. -I$(CMSSW_BASE)/src -std=c++0x
#CXXFLAGS  += -I. -I$(CMSSW_BASE)/src -std=c++0x -I$(LHAPDF_DATA_PATH)/../../include

# Optimization flag
CXXFLAGS += -g #-O3

# Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

# Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

INCLUDES += $(shell $(PYTHONCFG) --includes)

LIBS           = $(shell root-config --glibs)
LIBSTOPTAGGER += -L$(TTTDIR) -lTopTagger
#LHAPDFLIB      = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF

#OBJS       = $(patsubst %, $(ODIR)/%, $(OBJ))


PROGRAMS = MyAnalysis 


all: mkobj sampPyWrap $(PROGRAMS)

mkobj:
	@mkdir -p obj

#code to compile shared library to link samples to python
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SDIR)/samples.cc $(IDIR)/samples.h 
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.c
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

$(ODIR)/%.o : $(WDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

$(ODIR)/%.o : $(WDIR)/%.c
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

$(ODIR)/%.o : $(WDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

$(ODIR)/%.o : $(WDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IDIR) -I$(SDIR) -o $@ -c $<

MyAnalysis: $(ODIR)/MyAnalysis.o $(ODIR)/ExploreBackground.o $(ODIR)/ExploreTopTagger.o $(ODIR)/ExploreEventSelection.o $(ODIR)/NtupleClass.o $(ODIR)/Utility.o $(ODIR)/samples.o $(ODIR)/EventShapeVariables.o $(ODIR)/bdt_350to650_fwm10_jmtev_top6.o $(ODIR)/SetUpTopTagger.o
	$(LD) $^ $(LIBSTOPTAGGER) $(LIBS) -o $@

clean:
	rm -rf $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core $(ODIR)

-include $(ODIR)/*.d