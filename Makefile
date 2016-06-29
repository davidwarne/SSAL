#!/bin/make

CC = gcc
#OPTS = -O2 -std=gnu99 -D__SERIAL__ 
OPTS = -O2 -std=gnu99 -D__SERIAL__ -D__FLOAT64__ 
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

EXESRC = $(EXPDIR)/TestABC.c $(EXPDIR)/TestSSAL.c $(EXPDIR)/TestImportLSBML.c
SRC = $(IODIR)/cJSON.c $(SRCDIR)/SSAL.c $(ESSADIR)/segils.c $(ASSADIR)/satauls.c $(ASSADIR)/samlmcbs.c $(ASSADIR)/sactauls.c $(UTILDIR)/sumlnls.c  $(UTILDIR)/suhzds.c $(UTILDIR)/suarngs.c $(UTILDIR)/surngus.c $(UTILDIR)/surngexps.c $(UTILDIR)/surngpmfs.c $(UTILDIR)/surngpois.c  $(ESSADIR)/degils.c $(ASSADIR)/datauls.c $(ASSADIR)/damlmcbs.c $(ASSADIR)/dactauls.c $(UTILDIR)/dumlnls.c  $(UTILDIR)/duhzds.c $(UTILDIR)/durngus.c $(UTILDIR)/durngexps.c $(UTILDIR)/durngpmfs.c $(UTILDIR)/durngpois.c 
OBJS = $(SRC:.c=.o)
EXEOBJS=$(EXESRC:.c=.o)
EXE = $(EXESRC:.c=)
INC = -I ./include/  
LIBS = -lm -lgsl -lgslcblas

LIBDIR = lib
STATIC = $(LIBDIR)/lib$(NAME).a
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
