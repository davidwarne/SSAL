#!/bin/make

CC = gcc
#OPTS = -O2 -D__SERIAL__ 
OPTS = -g -pg  -D__SERIAL__ 
#OPTS = -g
SRCDIR=src
SRC =  $(SRCDIR)/TestSSAL.c $(SRCDIR)/SSAL.c $(SRCDIR)/segils.c $(SRCDIR)/satauls.c $(SRCDIR)/suhzds.c $(SRCDIR)/suarngs.c $(SRCDIR)/surngus.c $(SRCDIR)/surngexps.c $(SRCDIR)/surngpmfs.c $(SRCDIR)/surngpois.c 
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
