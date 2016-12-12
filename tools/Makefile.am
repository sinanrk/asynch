bin_PROGRAMS = asynch_graph asynch_peano

if HAVE_POSTGRESQL

bin_PROGRAMS += asynch_graph asynch_peano

asynch_createtables_SOURCES = createtables.c
asynch_createtables_LDADD = $(POSTGRESQL_LIBS)

asynch_deletetables_SOURCES = deletetables.c
asynch_deletetables_LDADD = $(POSTGRESQL_LIBS)

endif

asynch_graph_SOURCES = build_graph.c

asynch_peano_SOURCES = build_peano.c

AM_LDFLAGS = $(POSTGRESQL_LDFLAGS)
AM_CFLAGS = -I$(srcdir)/../src $(POSTGRESQL_CPPFLAGS)

endif