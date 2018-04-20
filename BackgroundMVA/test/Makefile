WDIR    = .
ODIR    = obj
IFWDIR  = $(CMSSW_BASE)/src/Framework/Framework/include
SFWDIR  = $(CMSSW_BASE)/src/Framework/Framework/src
IADIR   = $(CMSSW_BASE)/src/Analyzer/Analyzer/include
SADIR   = $(CMSSW_BASE)/src/Analyzer/Analyzer/src
SSDIR   = $(CMSSW_BASE)/src/SusyAnaTools/Tools
ITTDir  = $(CMSSW_BASE)/CfgParser/include
STTDir  = $(CMSSW_BASE)/CfgParser/src
TDIR    = $(CMSSW_BASE)/src

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


PROGRAMS = make_mva_training_tree_example tmva_train_example


all: mkobj sampPyWrap $(PROGRAMS)

mkobj:
	@mkdir -p obj

#code to compile shared library to link samples to python
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SFWDIR)/samples.cc $(IFWDIR)/samples.h 
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SFWDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

#Framework Repo
$(ODIR)/%.o : $(SFWDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SFWDIR)/%.c
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SFWDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SFWDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

#Analyzer Repo
$(ODIR)/%.o : $(SADIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SADIR)/%.c
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SADIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SADIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

#SusyAnaTools Repo
$(ODIR)/%.o : $(SSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SSDIR)/%.c
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(SSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

#TopTagger Repo
$(ODIR)/%.o : $(STTDir)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(STTDir)/%.c
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(STTDir)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(STTDir)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

#Working Directory 
$(ODIR)/%.o : $(WDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(WDIR)/%.c
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(WDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

$(ODIR)/%.o : $(WDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS)  -I$(TDIR) -I$(IFWDIR) -I$(SFWDIR) -I$(IADIR) -I$(SADIR) -I$(SSDIR) -I$(ITTDir) -I$(STTDir) -o $@ -c $<

#MyAnalysis: $(ODIR)/MyAnalysis.o $(ODIR)/AnalyzeBackground.o $(ODIR)/AnalyzeTopTagger.o $(ODIR)/AnalyzeEventSelection.o $(ODIR)/AnalyzeEventShape.o $(ODIR)/Analyze0Lep.o $(ODIR)/NTupleReader.o $(ODIR)/Utility.o $(ODIR)/samples.o $(ODIR)/EventShapeVariables.o $(ODIR)/bdt_350to650_fwm10_jmtev_top6.o $(ODIR)/SetUpTopTagger.o $(ODIR)/SATException.o
#	$(LD) $^ $(LIBSTOPTAGGER) $(LIBS) -o $@
#
#plot: $(ODIR)/plot.o
#	$(LD) $^ $(LIBS) -o $@

make_mva_training_tree_example: $(ODIR)/make_mva_training_tree_example.o
	$(LD) $^ $(LIBS) -o $@

tmva_train_example: $(ODIR)/tmva_train_example.o
	$(LD) $^ $(LIBS) -o $@

clean:
	rm -rf $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core $(ODIR)

-include $(ODIR)/*.d
