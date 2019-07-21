CC=gcc
head=./head/
CFLAG=-g -I$(head) 
sources=inc.o main.o update.o GridInit.o
a.out: $(sources)
	$(CC) -o a.out $(CFLAG) $(sources) -lm
update.o: update.c $(head)update.h
	$(CC) -c $(CFLAG) update.c
GridInit.o: GridInit.c $(head)GridInit.h
	$(CC) -c $(CFLAG) GridInit.c
inc.o: inc.c $(head)inc.h
	$(CC) -c $(CFLAG) inc.c 
main.o: main.c  $(head)fdtd.h
	$(CC) -c $(CFLAG) main.c

