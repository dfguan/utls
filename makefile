CC      =  gcc
CFLAGS  =  -g -Wall   
LDFLAGS = -lz

PROG = seq2sam  

.SUFFIXS:.c .o

all:$(PROG)

seq2sam: seq2sam.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	


