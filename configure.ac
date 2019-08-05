AC_INIT(cartogram, 1.0)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(cartogram_generator/main.c)

AM_INIT_AUTOMAKE([foreign])

dnl Get rid of the default -g -O2 CFLAGS used for gcc, since we want to use -O3.
save_CFLAGS=$CFLAGS
CFLAGS="-O3"
AC_PROG_CC
CFLAGS=$save_CFLAGS

AC_ARG_ENABLE([static], AC_HELP_STRING([--enable-static], [build a static binary (default is no)]), [
    CFLAGS="$CFLAGS -static"
    PKG_CONFIG_OPTS="--static"
], [PKG_CONFIG_OPTS=""]
)

AC_OPENMP

save_LIBS=$LIBS
save_CFLAGS=$CFLAGS
LIBS="$LIBS $OPENMP_CFLAGS"
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
AC_MSG_CHECKING([whether OpenMP will work in a package])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <omp.h>], [return omp_get_num_threads ();])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no]); AC_MSG_ERROR([A working OpenMP is required])]
)
LIBS=$save_LIBS
CFLAGS=$save_CFLAGS

dnl Checking for -lm
AC_CHECK_LIB(m,main, [], [AC_ERROR([A working math library is required])])

dnl Checking for -lfftw3
dnl AC_CHECK_HEADERS(fftw3.h, [], [AC_ERROR([A working fftw3 is required])])
dnl AC_CHECK_LIB(fftw3, fftw_plan_dft_r2c_2d, [], [AC_ERROR([A working fftw3 is required])])

AC_PATH_TOOL(PKG_CONFIG, pkg-config)

AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)
AC_MSG_CHECKING([for FFTW3])
if $PKG_CONFIG --exists fftw3; then
    FFTW3_CFLAGS=`$PKG_CONFIG $PKG_CONFIG_OPTS --cflags fftw3`
    FFTW3_LIBS=`$PKG_CONFIG $PKG_CONFIG_OPTS --libs fftw3`
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
fi

save_LIBS=$LIBS
save_CFLAGS=$CFLAGS
LIBS="$LIBS $FFTW3_LIBS"
CFLAGS="$CFLAGS $FFTW3_CFLAGS"
AC_MSG_CHECKING([whether FFTW3 will work in a package])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <fftw3.h>], [fftw_complex *in, *out;fftw_plan p;in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 50);out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 50);p = fftw_plan_dft_1d(50, in, out, FFTW_FORWARD, FFTW_ESTIMATE);fftw_destroy_plan(p);fftw_free(in); fftw_free(out);return 0;])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no]); AC_MSG_ERROR([A working FFTW3 is required])]
)
LIBS=$save_LIBS
CFLAGS=$save_CFLAGS

AC_SUBST(CJSON_CFLAGS)
AC_SUBST(CJSON_LIBS)
AC_MSG_CHECKING([for cJSON])
if $PKG_CONFIG --exists libcjson; then
    CJSON_CFLAGS=`$PKG_CONFIG $PKG_CONFIG_OPTS --cflags libcjson`
    CJSON_LIBS=`$PKG_CONFIG $PKG_CONFIG_OPTS --libs libcjson`
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
fi

save_LIBS=$LIBS
save_CFLAGS=$CFLAGS
LIBS="$LIBS $CJSON_LIBS"
CFLAGS="$CFLAGS $CJSON_CFLAGS"
AC_MSG_CHECKING([whether cJSON will work in a package])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <cjson/cJSON.h>], [cJSON *json = cJSON_Parse("");])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no]); AC_MSG_ERROR([A working cJSON is required])]
)
LIBS=$save_LIBS
CFLAGS=$save_CFLAGS

AC_OUTPUT(Makefile cartogram_generator/Makefile)