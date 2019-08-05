CC=clang 
head=./head/
CFLAG=  -I$(head) 
sources=inc.o main.o update.o GridInit.o abc.o update_pml.o
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
update_pml.o: update_pml.c $(head)update_pml.h
	$(CC) -c $(CFLAG) update_pml.c 
main.o: main.c  $(head)fdtd.h
	$(CC) -c $(CFLAG) main.c

