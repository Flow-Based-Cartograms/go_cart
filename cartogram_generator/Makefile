CC = gcc-7
CFLAGS = -O3 -fopenmp -Wall
DEPS = cartogram.h
LIBS = -lfftw3 -lm
OBJ = main.o cartogram.o integrate.o read_gen.o ps_figure.o fill_with_density.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

cartogram: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm *.o
