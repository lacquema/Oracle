#! /bin/bash

# Makefile to build libraries and executables

# Depedencies
COMPILF = gfortran
PYTHON3 = python3

# Others parameters
PAR = NO
UNIV = NO
ADD_FLAGS =
LIB_FLAGS = -O3 -c
ALG_FLAGS = -O3
DIR = .

# Induced directories
CODE_DIR = $(DIR)/Code

LIB_DIR = $(CODE_DIR)/lib
BIN_DIR = $(CODE_DIR)/bin
MCMC_DIR = $(CODE_DIR)/mcmc

MAIN_DIR = $(MCMC_DIR)/main
SUB_DIR = $(MCMC_DIR)/sub

# lib option 
LIB=mcmc
ifeq ($(UNIV),YES)
	LIB := $(LIB)_univ
endif
ifeq ($(PAR),YES)
	LIB_FLAGS+= -fopenmp
	ALG_FLAGS+= -fopenmp
	LIB := $(LIB)_par
endif

# Utilities
library_for:
	test ! -f $(LIB_DIR)/lib$(LIB).a || rm $(LIB_DIR)/lib$(LIB).a
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/utils.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mcmc.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/mrqfit.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/io.f 
ifeq ($(UNIV),YES)
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/kepu.f
else
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SUB_DIR)/kepellip.f
endif	
	ar -rcsv $(LIB_DIR)/lib$(LIB).a *.o
	rm *.o

library_py:
	$(PYTHON3) -m pip install -r $(DIR)/requirements.txt

# Algorithms avec ou sans _par
# astrom_mcmco:
# 	test ! -f $(BIN_DIR)/$@ || rm $(BIN_DIR)/$@
# 	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(MAIN_DIR)/$@.f -L$(LIB_DIR) -l$(LIB) -o $(BIN_DIR)/$@_par

# astrom_mcmco_univ:
# 	test ! -f $(BIN_DIR)/$@ || rm $(BIN_DIR)/$@
# 	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(MAIN_DIR)/$@.f -L$(LIB_DIR) -l$(LIB) -o $(BIN_DIR)/$@_par

astrom_%:
ifeq ($(PAR),NO)
	test ! -f $(BIN_DIR)/$@ || rm $(BIN_DIR)/$@
	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(MAIN_DIR)/$@.f -L$(LIB_DIR) -l$(LIB) -o $(BIN_DIR)/$@
else
	test ! -f $(BIN_DIR)/$@_par || rm $(BIN_DIR)/$@_par
	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(MAIN_DIR)/$@.f -L$(LIB_DIR) -l$(LIB) -o $(BIN_DIR)/$@_par
endif

compile:
	make library_for
	make library_for UNIV=YES
	make library_for PAR=YES
	make library_for UNIV=YES PAR=YES
	make astrom_mcmco
	make astrom_mcmco PAR=YES
	make astrom_univ_mcmco UNIV=YES
	make astrom_univ_mcmco UNIV=YES PAR=YES

all: 
	make compile
	make library_py

clean:
	rm *.o

cleanall: 
	clean
	rm *.a


