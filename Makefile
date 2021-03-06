#!/bin/make

CC = icc
#OPTS = -O2 -D__SERIAL__ -D__FLOAT64__
OPTS = -O2 -Wall -mkl=sequential -D__SERIAL__ -D__MKL__ 
#OPTS = -g  -fp-trap=common -Wall -mkl=sequential -D__SERIAL__ -D__MKL__ 
#OPTS = -pg -g -std=gnu99 -D__SERIAL__ 
#OPTS = -pg -g -std=gnu99 -D__SERIAL__ -D__FLOAT64__ 

NAME = ssal
#OPTS = -g   -D__SERIAL__ 
#OPTS = -g
EXPDIR=examples
SRCDIR=src
IODIR=$(SRCDIR)/io
UTILDIR=$(SRCDIR)/sequential/util
ESSADIR=$(SRCDIR)/sequential/essa
ASSADIR=$(SRCDIR)/sequential/assa
ODEDIR=$(SRCDIR)/sequential/ode
PDEDIR=$(SRCDIR)/sequential/pde

EXESRC = $(EXPDIR)/TestCRN.c $(EXPDIR)/TestCRN_cor.c $(EXPDIR)/TestSDE.c $(EXPDIR)/TestImportLSBML.c $(EXPDIR)/TestRNG.c $(EXPDIR)/TestODE.c $(EXPDIR)/TestPDE.c $(EXPDIR)/TestPDE_FKPP.c $(EXPDIR)/TestPDE_PF.c $(EXPDIR)/TestPhaseTypeDist.c $(EXPDIR)/TestDRW.c

SRC = $(IODIR)/cJSON.c $(SRCDIR)/SSAL.c $(UTILDIR)/suarngs.c $(ESSADIR)/degils.c $(ESSADIR)/demnrms.c $(ASSADIR)/datauls.c $(ASSADIR)/dactauls.c $(ASSADIR)/daectauls.c $(ASSADIR)/dalrws.c $(ASSADIR)/daems.c $(ASSADIR)/dacems.c $(UTILDIR)/duhzds.c $(UTILDIR)/durngus.c $(UTILDIR)/durngexps.c $(UTILDIR)/durngpmfs.c $(UTILDIR)/durngpois.c $(UTILDIR)/durngns.c $(UTILDIR)/durngmvns.c $(ODEDIR)/drk4s.c $(PDEDIR)/dbtcsfps.c
OBJS = $(SRC:.c=.o)
EXEOBJS=$(EXESRC:.c=.o)
EXE = $(EXESRC:.c=)
INC = -I ./include/  
LIBS = -lm

LIBDIR = lib
STATIC = $(LIBDIR)/lib$(NAME).a
#STATIC = $(LIBDIR)/lib$(NAME)_mkl.a
AR = ar
AROPTS = -rcvs
#PROFILE = -g -pg


.SUFFIXES: .c .o

all: $(EXE) 
	@echo Binaries $(EXE) created!

.c.o: 
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(STATIC): $(OBJS)
	$(AR) $(AROPTS) $(STATIC) $(OBJS)
#$(BIN): $(OBJS)
#	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) -lm 
#	@echo Binary created!!

$(EXE): $(EXEOBJS) $(STATIC) 
	$(CC) $(OPTS) -o $@  $(INC) $@.o  $(STATIC)  $(LIBS)
clean:
	set nonomatch; rm -f $(EXE) $(EXEOBJS) $(OBJS) $(STATIC)
