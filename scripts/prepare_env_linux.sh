# Going into source code directory
cd cartogram_generator

# Installing dependencies and generating Makefile
sudo apt-get update
sudo apt-get install libfftw3-3 libfftw3-dev && \
sudo apt-get install build-essential && \
echo "CC = gcc
CFLAGS = -O3 -fopenmp -Wall
DEPS = cartogram.h
LIBS = -lfftw3 -lm
OBJ = main.o cartogram.o ffb_integrate.o diff_integrate.o read_gen.o ps_figure.o fill_with_density.o

%.o: %.c \$(DEPS)
	\$(CC) -c -o \$@ \$< \$(CFLAGS)

cartogram: \$(OBJ)
	\$(CC) -o \$@ \$^ \$(CFLAGS) \$(LIBS)

.PHONY: clean

clean:
	rm *.o" > Makefile

# Returning to root directory
cd ..
