#! /bin/bash

# Makefile to build libraries and executables

# Parameters to input
COMPILF = /opt/homebrew/bin/gfortran

# Others parameters
PYTHON3 = /Users/lacquema/Oracle.env/bin/python3
PARALLEL = NO
ADD_FLAGS =
LIB_FLAGS = -O3 -c
ALG_FLAGS = -O3
DIR = .

# Induced directories
CODE_DIR = $(DIR)/Code

LIB_DIR = $(CODE_DIR)/lib
BIN_DIR = $(CODE_DIR)/bin
MCMC_DIR = $(CODE_DIR)/mcmc

ALG_DIR = $(MCMC_DIR)/main
SUB_DIR = $(MCMC_DIR)/sub

# Parallelization option 
ifeq ($(PARALLEL),NO)
LIB = $(LIB_DIR)/lib_mcmc.a
else
LIB_FLAGS += -fopenmp
ALG_FLAGS += -fopenmp
LIB = $(LIB_DIR)/lib_mcmc_par.a
endif

# Utilities
library_for:
	test ! -f $(LIB) || rm $(LIB)
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/utils.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mcmc.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mrqfit.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/io.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/kepellip.f
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/kepu.f
	ar -rv $(LIB) *.o
	rm *.o

library_py:
	$(PYTHON3) -m pip install -r $(DIR)/requirements.txt


# Algorithms
astrom_%:
ifeq ($(PARALLEL),NO)
	test ! -f $(BIN_DIR)/$@ || rm $(BIN_DIR)/$@
	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@
else
	test ! -f $(BIN_DIR)/$@_par || rm $(BIN_DIR)/$@_par
	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@_par
endif

compile:
	make library_for
	make astrom_mcmco
	make astrom_univ_mcmco
	make astrom_mcmco PARALLEL=YES
	make astrom_univ_mcmco PARALLEL=YES

all: 
	make compile
	make library_py

clean:
	rm *.o

cleanall: 
	clean
	rm *.a


