CC=gcc
CFLAGS= -fPIC -O3 -lm

all: fmftPy.o fmft.o nrutil.o
	$(CC) $(CFLAGS) -shared -o libfmft.so $^

samfmft: samfmft.c nrutil.o fmft.o
	gcc -o $@ $^ -lm
main_fmft: main_fmft.c nrutil.o fmft.o
	gcc -o $@ $^  -lm
%.o: %.c 
	gcc -fPIC -c $< 
clean:
	rm -f *.o *.mod 
