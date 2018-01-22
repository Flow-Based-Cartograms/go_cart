# Going into source code directory
cd cartogram_generator

# Install dependencies and generate Makefile
brew install fftw && \
gcc_version=$(brew install gcc 2>&1 | tee /dev/tty | egrep -o "(gcc)(\D*)(\d+)" | head -1 | sed -E 's|^[^0-9]*||') && \
echo "CC = gcc-$gcc_version
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

# Linking fftw
brew link fftw | echo "Linking FFTW"

# Returning to root directory
cd ..
