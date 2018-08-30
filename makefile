# Copyright (c) 2014-2018 Michael C. Heiber
# This source file is part of the Ising_OPV project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

COMPILER = mpicxx
FLAGS = -Wall -Wextra -O3 -std=c++11
OBJS = src/main.o src/Morphology.o src/Lattice.o src/Utils.o tinyxml2/tinyxml2.o

Ising_OPV.exe : $(OBJS)
	$(COMPILER) $(FLAGS) $(OBJS) -o Ising_OPV.exe

src/main.o : src/main.cpp src/Morphology.h src/Lattice.h src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@
	
src/Morphology.o : src/Morphology.h src/Morphology.cpp src/Lattice.h src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@

src/Lattice.o : src/Lattice.h src/Lattice.cpp src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@
	
src/Utils.o : src/Utils.h src/Utils.cpp
	$(COMPILER) $(FLAGS) -c $< -o $@
	
tinyxml2/tinyxml2.o : tinyxml2/tinyxml2.h tinyxml2/tinyxml2.cpp
	$(COMPILER) $(FLAGS) -c $< -o $@

clean:
	\rm src/*.o tinyxml/*.o *~ Ising_OPV.exe