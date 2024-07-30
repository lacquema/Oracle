#! /bin/sh
#
# Makefile to build simulations
#

SHELL = /bin/sh
ARCH = $(shell uname)
DEBUGGING=NO
PARALLEL=NO
DEBUG = -g
FFLAGS = -O3
LDFLAGS = -O3 
#ulimit -s unlimited
#GAG_EXEC_DIR = /home/ciment/lacquema/gildas-exe-feb09c/x86_64-redhatServer-ifort
#GAG_PATH = /home/ciment/lacquema/gildas-exe-feb09c/etc
#LIBGAG = -L$(GAG_EXEC_DIR)/lib -lgio -lginc
#LOCAL_GAG = /home/ciment/lacquema/greg

BIN_DIR = ./bin

ifeq ($(DEBUGGING),YES)
FFLAGS += $(DEBUG)
LDFLAGS = $(DEBUG)
endif
ifeq ($(PARALLEL),YES)
FFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif
FC = gfortran

%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

% : %.o
ifeq ($(PARALLEL),YES)
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$(join $@,p) $?
else
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$@ $?
endif
	rm $?

clean:
	rm -f *.o *~ core

# Simulation

MCMC = mcmc/utils.o mcmc/mcmc.o mcmc/mrqfit.o mcmc/radvel.o mcmc/astrometry.o
astrom_mcmco: $(MCMC) astrom_mcmco.o



