CC=gcc
CFLAG=-lm 
a.out: main.c
	$(CC) $(CFLAG) main.c
