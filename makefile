CC      =  gcc
CFLAGS  =  -g -Wall  -D TEST 
LDFLAGS = -lz

PROG = seq2sam eutls 10x_trim

.SUFFIXS:.c .o

all:$(PROG)

seq2sam: seq2sam.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

eutls: eutls.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

10x_trim: bseq.o 10x_trim.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	


