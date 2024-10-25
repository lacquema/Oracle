#! /bin/bash

# Makefile to build libraries and executables

# Parameters
COMPILF = /opt/homebrew/bin/gfortran
PYTHON3 = /Users/lacquema/Oracle.env/bin/python3
PARALLEL = NO
ADD_FLAGS =
LIB_FLAGS = -O3 -c
ALG_FLAGS = -O3
DIR = .

# Induced directories
CODE_DIR = $(DIR)/code

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
library:
	test ! -f $(LIB) || rm $(LIB)
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/utils.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mcmc.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mrqfit.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/io.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/kepellip.f
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/kepu.f
	ar -rv $(LIB) *.o
	rm *.o

packages:
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

all: 
	make packages
	make library
	make astrom_mcmco
	make astrom_univ_mcmco

clean:
	rm *.o

cleanall: 
	clean
	rm *.a


