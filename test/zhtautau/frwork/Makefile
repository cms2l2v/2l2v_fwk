CXX=g++
CXXFLAGS=-O -Wall -fPIC
RINT=-DWITHRINT
THREAD_LIB=-lpthread

#use TRint or TApplication
CXXFLAGS += $(RINT)

## Libs
ROOTCXXFLAGS = $(shell root-config --cflags)
ROOTLDFLAGS  = $(shell root-config --ldflags)
ROOTLIBS     = $(shell root-config --libs)
ROOTGLIBS    = $(shell root-config --glibs)

LIBS  = ${THREAD_LIB} -lcurses 
LIBS += $(ROOTLIBS)

CXXFLAGS += $(ROOTCXXFLAGS)

#objects
PLOTT_OBJ  = Analyzer.o
PLMAK_OBJ  = FRAnalyzer.o
##Build
all: FRAnalyzer Analyzer

Analyzer:	$(PLOTT_OBJ)
		@echo ""
		@echo "Compiling objects ($(PLOTT_OBJ)) ... "
		@$(CXX) $(PLOTT_OBJ) $(ROOTLDFLAGS) -o Analyzer FRAnalyzer.o -g $(LIBS) -lm -lEG

FRAnalyzer:	$(PLMAK_OBJ)
		@echo ""
		@echo "Compiling objects ($(PLMAK_OBJ)) ... "
		@$(CXX) -c $(PLMAK_OBJ) $(ROOTLDFLAGS) -o -g $(LIBS) -lm -lEG

## Clean
clean:
	@echo "Cleaning path ... "
	rm -rf *.o
	rm -rf Analyzer
	rm -rf *~
