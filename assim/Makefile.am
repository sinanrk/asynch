noinst_LIBRARIES = libforecaster.a libassim.a

libforecaster_a_SOURCES = forecaster_methods.c
libassim_a_SOURCES = assim_ls_methods.c assim_models.c

bin_PROGRAMS = assim_least_square_linear #assimpersis assimpersis_end

include_HEADERS = forecaster_methods.h

assim_least_square_linear_SOURCES = assim_least_squares_linear.c
assim_least_square_linear_LDADD = libforecaster.a libassim.a $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS) $(METIS_LIBS) $(PETSC_LIBS)

#assimpersis_SOURCES = assimpersis .c
#assimpersis_LDADD = libforecaster.a libassim.a $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS)

#assimpersis_end_SOURCES = assimpersis_end.c
#assimpersis_end_LDADD = libforecaster.a libassim.a $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS) $(HDF5_LIBS) $(SSH2_LIBS)

AM_LDFLAGS = $(POSTGRESQL_LDFLAGS) $(HDF5_LDFLAGS) $(SSH2_LDFLAGS)
AM_CFLAGS = -I$(srcdir)/../src $(POSTGRESQL_CPPFLAGS) $(HDF5_CFLAGS) $(SSH2_CPPFLAGS) $(METIS_CFLAGS) $(PETSC_CFLAGS)
