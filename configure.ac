AC_INIT([asynch], [1.3.2], [samuel-debionne@uiowa.edu])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Check for programms
AX_PROG_CC_MPI([], [], [AC_MSG_FAILURE([MPI compiler requested, but couldn't find MPI.])] )
AC_PROG_CC_STDC #C99
AM_PROG_CC_C_O  #automake < 1.14
AM_PROG_AR
AC_PROG_RANLIB

# Checks for libraries
AX_CHECK_ZLIB
AX_LIB_HDF5([serial])
AX_LIB_POSTGRESQL
AX_BLAS
AX_LAPACK

# Checks for packages
PKG_CHECK_MODULES([CHECK], [check >= 0.9.8], [AC_DEFINE([HAVE_LIBCHECK], [1], [Defined if you have libcheck support] )],[HAVE_LIBCHECK=0])
PKG_CHECK_MODULES([METIS], [metis >= 5.0.0], [AC_DEFINE([HAVE_METIS], [1], [Defined if you have METIS support] )], [HAVE_METIS=0])
PKG_CHECK_MODULES([PETSC], [PETSc >= 3.7.0], [AC_DEFINE([HAVE_PETSC], [1], [Defined if you have PETSc support] )], [HAVE_PETSC=0])

AM_CONDITIONAL([HAVE_POSTGRESQL], [test "$HAVE_POSTGRESQL" -eq 1])
AM_CONDITIONAL([HAVE_PETSC], [test "$HAVE_PETSC" -eq 1])

# Generate config header
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile forecasters/Makefile tests/Makefile tools/Makefile])

AC_OUTPUT
