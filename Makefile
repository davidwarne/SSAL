#!/bin/make

CC = gcc
OPTS = -O2 -std=gnu99 -D__SERIAL__ 
#OPTS = -g   -D__SERIAL__ 
#OPTS = -g
EXPDIR=examples
SRCDIR=src
UTILDIR=$(SRCDIR)/sequential/util
ESSADIR=$(SRCDIR)/sequential/essa
ASSADIR=$(SRCDIR)/sequential/assa

EXESRC = $(EXPDIR)/TestABC.c $(EXPDIR)/TestSSAL.c
SRC = $(SRCDIR)/SSAL.c $(ESSADIR)/segils.c $(ASSADIR)/satauls.c $(ASSADIR)/samlmcbs.c $(ASSADIR)/sactauls.c $(UTILDIR)/sumlnls.c  $(UTILDIR)/suhzds.c $(UTILDIR)/suarngs.c $(UTILDIR)/surngus.c $(UTILDIR)/surngexps.c $(UTILDIR)/surngpmfs.c $(UTILDIR)/surngpois.c 
OBJS = $(SRC:.c=.o)
EXEOBJS=$(EXESRC:.c=.o)
EXE = $(EXESRC:.c=)
INC = -I ./include/ 
BIN = TestABC
LIBS = -lm
#PROFILE = -g -pg


.SUFFIXES: .c .o

all: $(EXE) 
	@echo Binaries $(EXE) created!

.c.o: 
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

#$(BIN): $(OBJS)
#	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) -lm 
#	@echo Binary created!!

$(EXE): $(OBJS) $(EXEOBJS)
	$(CC) $(OPTS) -o $@  $(INC) $@.o $(OBJS) $(LIBS)
clean:
	set nonomatch; rm -f $(EXE) $(EXEOBJS) $(OBJS)
