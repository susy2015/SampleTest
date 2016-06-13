IDIR = .
ODIR = obj
SDIR = ../../SusyAnaTools/Tools/

CXX = g++
CXXFLAGS += -I. -std=c++0x -I${CMSSW_BASE}/src/ -I$(LHAPDF_DATA_PATH)/../../include
## Optimization flag
CXXFLAGS += -g
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXDEPFLAGS = -MMD -MP
LD = g++
LDFLAGS =
LIBS = $(shell root-config --glibs)

LHAPDFLIB = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF
#OBJS = $(patsubst %, $(ODIR)/%, $(OBJ))

PROGRAMS = MC Data CS HadTauMC LLMC PrintSearchBins

all: mkobj sampPyWrap $(PROGRAMS)

mkobj:	
	@mkdir -p obj

#code to compile shared library to link samples to python
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SDIR)/samples.cc $(SDIR)/samples.h
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<


MC: $(ODIR)/MC.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

HadTauMC: $(ODIR)/HadTauMC.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

LLMC: $(ODIR)/LLMC.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

CS: $(ODIR)/CS.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

Data: $(ODIR)/Data.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

PrintSearchBins: $(ODIR)/PrintSearchBins.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(PROGRAMS) core

-include $(ODIR)/*.d