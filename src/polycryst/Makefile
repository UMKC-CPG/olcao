
CC = gcc
FL = -g -Wall
LB = -lm
OJ = -c
rev := $(shell date +%y%m%d)

all: polycryst

polycryst: polycryst.c functions.o
	${CC} ${FL} -o polycryst polycryst.c functions.o ${LB}

functions.o: functions.c
	${CC} ${FL} ${OJ} functions.c

clean:
	rm -f polycryst *.o

tidy:
	rm -f polycryst.in polycryst.dat euler_parameters *.o site siteq input

run: polycryst tidy
	./polycryst

backup:
	tar czvf ./codes/polycryst$(rev).tar.gz *.c *.h Makefile

