COMPILER = mpicxx
FLAGS = -Wall -Wextra -O3 -std=c++11
OBJS = main.o Morphology.o Lattice.o Utils.o

Ising_OPV.exe : $(OBJS)
	$(COMPILER) $(FLAGS) $(OBJS) -o Ising_OPV.exe

main.o : main.cpp Morphology.h Lattice.h Utils.h
	$(COMPILER) $(FLAGS) -c main.cpp
	
Morphology.o : Morphology.h Morphology.cpp Lattice.h Utils.h
	$(COMPILER) $(FLAGS) -c Morphology.cpp

Lattice.o : Lattice.h Lattice.cpp Utils.h
	$(COMPILER) $(FLAGS) -c Lattice.cpp
	
Utils.o : Utils.h Utils.cpp
	$(COMPILER) $(FLAGS) -c Utils.cpp
	
clean:
	\rm *.o *~ Ising_OPV.exe