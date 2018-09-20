# Copyright (c) 2014-2018 Michael C. Heiber
# This source file is part of the Ising_OPV project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

ifeq ($(lastword $(subst /, ,$(CXX))),g++)
	FLAGS += -Wall -Wextra -O3 -std=c++11 -I. -Isrc
endif
ifeq ($(lastword $(subst /, ,$(CXX))),pgc++)
	FLAGS += -O2 -fastsse -Mvect -std=c++11 -Mdalign -Munroll -Mipa=fast -Kieee -m64 -I. -Isrc
endif

COMPILER = mpicxx
OBJS = src/Morphology.o src/Lattice.o src/Parameters.o src/Utils.o tinyxml2/tinyxml2.o

all : Ising_OPV.exe
ifndef FLAGS
	$(error Valid compiler not detected.)
endif

Ising_OPV.exe : src/main.o $(OBJS)
	$(COMPILER) $(FLAGS) $^ -o $@

src/main.o : src/main.cpp src/Morphology.h src/Lattice.h src/Parameters.h src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@
	
src/Morphology.o : src/Morphology.cpp src/Morphology.h src/Lattice.h src/Parameters.h src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@

src/Lattice.o : src/Lattice.cpp src/Lattice.h src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@

src/Parameters.o : src/Parameters.cpp src/Parameters.h src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@	

src/Utils.o : src/Utils.cpp src/Utils.h
	$(COMPILER) $(FLAGS) -c $< -o $@
	
tinyxml2/tinyxml2.o : tinyxml2/tinyxml2.cpp tinyxml2/tinyxml2.h
	$(COMPILER) $(FLAGS) -c $< -o $@

#
# Testing Section using googletest
#

ifndef FLAGS
	$(error Valid compiler not detected.)
endif
GTEST_DIR = googletest/googletest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
ifeq ($(lastword $(subst /, ,$(CXX))),g++)
	GTEST_FLAGS = -isystem $(GTEST_DIR)/include -pthread
endif
ifeq ($(lastword $(subst /, ,$(CXX))),pgc++)
	GTEST_FLAGS = -I$(GTEST_DIR)/include
endif

test_coverage : FLAGS = -fprofile-arcs -ftest-coverage -std=c++11 -Wall -Wextra -I. -Isrc
test_coverage : test/Ising_OPV_tests.exe test/Ising_OPV_MPI_tests.exe

test : test/Ising_OPV_tests.exe test/Ising_OPV_MPI_tests.exe
	
test/Ising_OPV_tests.exe : test/test.o test/gtest-all.o $(OBJS)
	mpicxx $(GTEST_FLAGS) $(FLAGS) $^ -lpthread -o $@

test/gtest-all.o : $(GTEST_SRCS_)
	mpicxx $(GTEST_FLAGS) -I$(GTEST_DIR) $(FLAGS) -c $(GTEST_DIR)/src/gtest-all.cc -o $@
			
test/test.o : test/test.cpp $(GTEST_HEADERS) $(OBJS)
	mpicxx $(GTEST_FLAGS) $(FLAGS) -c $< -o $@
	
test/Ising_OPV_MPI_tests.exe : test/test_mpi.o test/gtest-all.o $(OBJS)
	mpicxx $(GTEST_FLAGS) $(FLAGS) -lpthread $^ -o $@

test/test_mpi.o : test/test_mpi.cpp $(GTEST_HEADERS) $(OBJS)
	mpicxx $(GTEST_FLAGS) $(FLAGS) -c $< -o $@
	
clean:
	\rm src/*.o tinyxml/*.o *~ Ising_OPV.exe src/*.gcno* src/*.gcda test/*.o test/*.gcno* test/*.gcda test/Ising_OPV_tests.exe test/Ising_OPV_MPI_tests.exe