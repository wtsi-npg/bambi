AM_CPPFLAGS = $(HTSLIB_CPPFLAGS) -std=gnu99
AM_LDFLAGS = -rdynamic $(HTSLIB_LDFLAGS)

bin_PROGRAMS = bambi 
bambi_SOURCES = bambi.c bclfile.c  decode.c   filterfile.c  hts_addendum.c  i2b.c  posfile.c array.c select.c \
                bambi.h bclfile.h             filterfile.h  hts_addendum.h         posfile.h array.h \
                parse.c bamit.c    chrsplit.c read2tags.c   \
                parse.h bamit.h
nobase_include_HEADERS = cram/cram_samtools.h cram/pooled_alloc.h cram/sam_header.h cram/string_alloc.h
bambi_CFLAGS = -I/usr/include/libxml2
bambi_LDADD = $(HTSLIB_HOME)/lib/libhts.a -ldl -lxml2 -lz -llzma -lbz2 -lpthread
