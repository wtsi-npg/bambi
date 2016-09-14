AM_CPPFLAGS = $(HTSLIB_CPPFLAGS)
AM_LDFLAGS = $(HTSLIB_LDFLAGS)

check_PROGRAMS = t_i2b
t_i2b_SOURCES = t_i2b.c ../../src/posfile.c ../../src/bclfile.c ../../src/filterfile.c ../../src/array.c
t_i2b_CFLAGS = -I$(top_srcdir)/src -I/usr/include/libxml2 -DDATA_DIR=$(abs_srcdir)/data
t_i2b_LDADD = $(HTSLIB_HOME)/lib/libhts.a -ldl -lxml2 -lz -llzma -lbz2 -lpthread

EXTRA_DIST = data
