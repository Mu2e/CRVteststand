CXX = g++
CXXFLAGS = -Wall
INCLIBFLAGS = `root-config --cflags --glibs` -lstdc++fs
#INCLIBFLAGS = `root-config --cflags --glibs` -lstdc++fs -I$(CLHEP_INC)

#TARGETS = parserCrv calibCrv calibLDMX recoCrv recoCrvPulseHeight recoCrvBeta recoCrvShortTail reflectedPulses mergeTestbeam merge tracklengths
#TARGETS = parserCrv calibCrv recoCrv reflectedPulses mergeTestbeam merge
TARGETS = parserCrv calibCrv recoCrv

all: $(TARGETS)

$(TARGETS): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCLIBFLAGS)

clean:
	rm $(TARGETS)
