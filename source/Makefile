# including other makefiles.

Homeinc = ../ARPACK

include $(Homeinc)/ARmake.inc


Targets = DVCI

all: DVCI Transitions FinalVCI

DVCI: DVCI.o Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o
	$(CPP) $(CPP_FLAGS) -o DVCI DVCI.o \
	Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o $(ALL_LIBS);
	echo "\n****DVCI compiled*****\n";

Transitions: Transitions.o Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o
	$(CPP) $(CPP_FLAGS) -o Transitions Transitions.o \
	Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o $(ALL_LIBS);
	echo "\n*****Transitions compiled*****\n";

FinalVCI: FinalVCI.o Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o
	$(CPP) $(CPP_FLAGS) -o FinalVCI FinalVCI.o \
	Graph.o ReadData.o Shared.o Assemblage.o Solver2.o Generator.o $(ALL_LIBS);
	echo "\n****FinalVCI compiled*****\n";

Graph.o: Graph.cc Graph.h Shared.o
	$(CPP) $(CPP_FLAGS) -c Graph.cc -o Graph.o

Generator.o: Generator.cc Generator.h Shared.o  Graph.o 
	$(CPP) $(CPP_FLAGS) -c Generator.cc -o Generator.o

ReadData.o: ReadData.h ReadData.cc Shared.o
	$(CPP) $(CPP_FLAGS) -c ReadData.cc -o ReadData.o

Assemblage.o: Assemblage.h Assemblage.cc Graph.o Shared.o
	$(CPP) $(CPP_FLAGS)  -c Assemblage.cc -o Assemblage.o

Shared.o: Shared.h Shared.cc
	$(CPP) $(CPP_FLAGS) -c Shared.cc -o Shared.o

Solver2.o: Solver2.h Solver2.cc Shared.o 
	$(CPP) $(CPP_FLAGS) -c Solver2.cc -o Solver2.o

DVCI.o: DVCI.cc Solver2.o Assemblage.o Basic.h Shared.o Graph.o Generator.o
	$(CPP) $(CPP_FLAGS) -c $< -o $@

FinalVCI.o: FinalVCI.cc Basic.h Shared.o ReadData.o Solver2.o
	$(CPP) $(CPP_FLAGS) -c $< -o $@

Transitions.o: Transitions.cc Solver2.o Assemblage.o Basic.h Shared.o Graph.o Generator.o
	$(CPP) $(CPP_FLAGS) -c $< -o $@

clean:
	touch $(Targets); rm -rf *.o;

mrproper: clean
	rm -rf main

