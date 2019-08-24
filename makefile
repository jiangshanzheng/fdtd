CC=clang 
head=./head/
CFLAG= -O3 -I$(head) 
sources=inc.o main.o update.o GridInit.o abc.o filter.o
a.out: $(sources)
	$(CC) -o a.out $(CFLAG) $(sources) -lm
update.o: update.c $(head)update.h
	$(CC) -c $(CFLAG) update.c
GridInit.o: GridInit.c $(head)GridInit.h
	$(CC) -c $(CFLAG) GridInit.c
inc.o: inc.c $(head)inc.h
	$(CC) -c $(CFLAG) inc.c 
abc.o: abc.c $(head)abc.h
	$(CC) -c $(CFLAG) abc.c 
filter.o: filter.c $(head)filter.h
	$(CC) -c $(CFLAG) filter.c 
main.o: main.c  $(head)fdtd.h
	$(CC) -c $(CFLAG) main.c
clean:
	rm *.o
