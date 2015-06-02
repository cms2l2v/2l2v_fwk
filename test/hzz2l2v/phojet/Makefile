#
#  Xin Shi <Xin.Shi@cern.ch>
#  2014.03.13 
#

CC = g++
GCCFLAGS  = -Wall -g 

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))


DRAWFLAGS=$(GCCFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -lHistPainter 


all: drawPt drawEff
	@echo "Build successful."

drawPt: drawPt.cc
	$(CC) $(DRAWFLAGS) $< -o $@

drawEff: drawEff.cc
	$(CC) $(DRAWFLAGS) $< -o $@

clean:
	rm -rf *~  *.dSYM 
	rm -f drawPt drawEff


