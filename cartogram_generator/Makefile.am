bin_PROGRAMS = cartogram
cartogram_SOURCES = cartogram.h cartogram.c diff_integrate.c ffb_integrate.c fill_with_density.c main.c ps_figure.c read_map.c process_json.c

cartogram_CFLAGS = -O3 -Wall $(OPENMP_CFLAGS) $(FFTW3_CFLAGS) $(CJSON_CFLAGS)
cartogram_LDADD = $(OPENMP_CFLAGS) $(FFTW3_LIBS) $(CJSON_LIBS)
