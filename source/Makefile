# including other makefiles.

Homeinc = ../ARPACK

include $(Homeinc)/ARmake.inc


Targets = DVCI

all: DVCI

DVCI: DVCI.o Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o 
	$(CPP) $(CPP_FLAGS)  -o DVCI DVCI.o \
	Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o $(ALL_LIBS) 

Graph.o: Graph.cc Graph.h Shared.cc
	$(CPP) $(CPP_FLAGS)  -c Graph.cc -o Graph.o

Generator.o: Generator.cc Shared.cc Generator.h Graph.cc 
	$(CPP) $(CPP_FLAGS) -c Generator.cc -o Generator.o

ReadData.o: ReadData.h Shared.cc ReadData.cc
	$(CPP) $(CPP_FLAGS) -c ReadData.cc -o ReadData.o

Assemblage.o: Assemblage.h Shared.cc Assemblage.cc Graph.cc
	$(CPP) $(CPP_FLAGS)  -c Assemblage.cc -o Assemblage.o

Shared.o: Shared.h Shared.cc
	$(CPP) $(CPP_FLAGS) -c Shared.cc -o Shared.o

Solver2.o: Solver2.h Solver2.cc Shared.cc 
	$(CPP) $(CPP_FLAGS) -c Solver2.cc -o Solver2.o

DVCI.o: DVCI.cc Solver2.cc Assemblage.cc Basic.h Shared.cc Graph.cc Generator.cc
	$(CPP) $(CPP_FLAGS)  -c $< -o $@

clean:
	touch $(Targets); rm -rf *.o;

mrproper: clean
	rm -rf main

