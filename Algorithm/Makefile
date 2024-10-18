#! /bin/bash

# Makefile to build libraries and executables

# Parameters
F = gfortran
PARALLEL = NO
ADD_FLAGS =
LIB_FLAGS = -O3 -c
ALG_FLAGS = -O3
DIR = .

# Induced directories
MCMC_DIR = $(DIR)/mcmc
LIB_DIR = $(DIR)/lib
BIN_DIR = $(DIR)/bin
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
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/utils.f 
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mcmc.f 
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mrqfit.f 
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/io.f 
	ar -rv $(LIB) *.o
	rm *.o

# Algorithms
astrom_mcmco:
ifeq ($(PARALLEL),NO)
	$(F) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@
else
	$(F) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@_par
endif

all: 
	make library
	make astrom_mcmco

clean:
	rm *.o

cleanall: 
	clean
	rm *.a


