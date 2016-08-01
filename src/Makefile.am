AM_CPPFLAGS = $(HTSLIB_CPPFLAGS)
AM_LDFLAGS = -rdynamic $(HTSLIB_LDFLAGS)

bin_PROGRAMS = bambi 
bambi_SOURCES = bambi.c bclfile.c  decode.c  filterfile.c  hts_addendum.c  i2b.c  posfile.c
bambi_CFLAGS = -I/usr/include/libxml2
bambi_LDADD = $(HTSLIB_HOME)/lib/libhts.a -ldl -lxml2 -lz -llzma -lbz2 -lpthread
