COMPILER = mpicxx
FLAGS = -Wall -Wextra -O3 -std=c++11
OBJS = main.o Morphology.o Lattice.o Utils.o tinyxml2.o

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
	
tinyxml2.o : tinyxml2/tinyxml2.h tinyxml2/tinyxml2.cpp
	$(COMPILER) $(FLAGS) -c tinyxml2/tinyxml2.cpp

clean:
	\rm *.o *~ Ising_OPV.exe