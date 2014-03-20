Ising_OPV.exe : main.cpp Morphology.cpp
	mpicxx -o Ising_OPV.exe main.cpp Morphology.cpp -Wall -Wextra -O2