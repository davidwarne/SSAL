#!/bin/make

CC = gcc
#OPTS = -O2 -D__SERIAL__ 
OPTS = -g -pg  -D__SERIAL__ 
#OPTS = -g 
SRC =  TestSSAL.c SSAL.c segils.c suhzds.c 
OBJS = $(SRC:.c=.o)
INC =  
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
