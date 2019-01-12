CC      =  gcc
CFLAGS  =  -g -Wall   
LDFLAGS = -lz

PROG = fast2sam  

.SUFFIXS:.c .o

all:$(PROG)

fast2sam: fast2sam.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	


