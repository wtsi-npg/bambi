AM_CPPFLAGS = $(HTSLIB_CPPFLAGS)
AM_LDFLAGS = $(HTSLIB_LDFLAGS)

bin_PROGRAMS = t_i2b
t_i2b_SOURCES = t_i2b.c ../../src/posfile.c ../../src/bclfile.c ../../src/filterfile.c
t_i2b_CFLAGS = -I$(top_srcdir)/src -I/usr/include/libxml2 -DDATA_DIR=$(abs_srcdir)
t_i2b_LDADD = $(HTSLIB_HOME)/lib/libhts.a -ldl -lxml2 -lz -llzma -lbz2 -lpthread
