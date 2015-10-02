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
SRC =  $(EXPDIR)/TestSSAL.c $(SRCDIR)/SSAL.c $(ESSADIR)/segils.c $(ASSADIR)/satauls.c $(UTILDIR)/suhzds.c $(UTILDIR)/suarngs.c $(UTILDIR)/surngus.c $(UTILDIR)/surngexps.c $(UTILDIR)/surngpmfs.c $(UTILDIR)/surngpois.c 
OBJS = $(SRC:.c=.o)
INC = -I ./include/ 
BIN = TestSSAL
LIBS = -lm
#PROFILE = -g -pg


.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) -lm 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) $(OBJS)
