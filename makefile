CC=gcc
CFLAG=-O3 -lm 
a.out: main.c
	$(CC) $(CFLAG) main.c
