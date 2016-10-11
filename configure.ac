AC_INIT([asynch], [1.0.1], [felipe-quintero@uiowa.edu])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([tar-ustar foreign -Wall -Werror])

# Check for programms
AX_PROG_CC_MPI([], [], [AC_MSG_FAILURE([MPI compiler requested, but couldn't find MPI.])] )
AC_PROG_CC_STDC #C99
AC_PROG_RANLIB

# Checks for libraries
AX_CHECK_ZLIB
AX_LIB_POSTGRESQL
AX_LIB_HDF5([serial])
AX_BLAS
AX_LAPACK

# Checks for packages
PKG_CHECK_MODULES([CHECK], [check >= 0.9.8], [AC_DEFINE([HAVE_LIBCHECK], [1], [Defined if you have libcheck support] )],[HAVE_LIBCHECK=0])
PKG_CHECK_MODULES([METIS], [check >= 5.0.0], [AC_DEFINE([HAVE_METIS], [1], [Defined if you have METIS support] )], [HAVE_METIS=0])
PKG_CHECK_MODULES([PETSC], [check >= 3.7.0], [AC_DEFINE([HAVE_PETSC], [1], [Defined if you have PETSC support] )], [HAVE_PETSC=0])

# Generate config header
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile assim/Makefile forecasters/Makefile tests/Makefile tools/Makefile])

AC_OUTPUT
