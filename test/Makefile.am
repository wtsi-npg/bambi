AM_CPPFLAGS = $(HTSLIB_CPPFLAGS) -std=gnu99 -Werror
AM_LDFLAGS = $(HTSLIB_LDFLAGS)

check_PROGRAMS = t_read2tags t_array t_bamit t_chrsplit t_select t_bclfile t_decode t_filterfile t_posfile t_i2b

t_read2tags_SOURCES = t_read2tags.c $(top_srcdir)/src/read2tags.c ../src/array.c ../src/bamit.c ../src/parse.c ../src/hts_addendum.c
t_read2tags_CFLAGS = -I$(top_srcdir)/src -DDATA_DIR=$(abs_srcdir)/data
t_read2tags_LDADD = $(HTSLIB_HOME)/lib/libhts.a -lz -ldl -lxml2 -lpthread -llzma -lbz2 -lcurl -lcrypto

t_select_SOURCES = t_select.c $(top_srcdir)/src/select.c ../src/array.c ../src/bamit.c ../src/parse.c ../src/hts_addendum.c
t_select_CFLAGS = -I$(top_srcdir)/src -DDATA_DIR=$(abs_srcdir)/data
t_select_LDADD = $(HTSLIB_HOME)/lib/libhts.a -lz -ldl -lxml2 -lpthread -llzma -lbz2 -lcurl -lcrypto

t_chrsplit_SOURCES = t_chrsplit.c $(top_srcdir)/src/chrsplit.c ../src/array.c ../src/bamit.c ../src/parse.c ../src/hts_addendum.c
t_chrsplit_CFLAGS = -I$(top_srcdir)/src -DDATA_DIR=$(abs_srcdir)/data
t_chrsplit_LDADD = $(HTSLIB_HOME)/lib/libhts.a -lz -ldl -lxml2 -lpthread -llzma -lbz2 -lcurl -lcrypto

t_bamit_SOURCES = t_bamit.c $(top_srcdir)/src/bamit.c 
t_bamit_CFLAGS = -I$(top_srcdir)/src -DDATA_DIR=$(abs_srcdir)/data
t_bamit_LDADD = $(HTSLIB_HOME)/lib/libhts.a -lz -ldl -lxml2 -lpthread -llzma -lbz2 -lcurl -lcrypto

t_array_SOURCES = t_array.c $(top_srcdir)/src/array.c 
t_array_CFLAGS = -I$(top_srcdir)/src -DDATA_DIR=$(abs_srcdir)/data
t_array_LDADD = $(HTSLIB_HOME)/lib/libhts.a -lz -ldl -lxml2 -lpthread -llzma -lbz2 -lcurl -lcrypto

t_bclfile_SOURCES = t_bclfile.c
t_bclfile_CFLAGS = -DDATA_DIR=$(abs_srcdir)/data
t_bclfile_LDADD = -lz -ldl -lxml2 -lpthread -lcurl -lcrypto

t_decode_SOURCES = t_decode.c ../src/array.c ../src/bamit.c
t_decode_CFLAGS = -I$(top_srcdir)/src -DDATA_DIR=$(abs_srcdir)/data
t_decode_LDADD = $(HTSLIB_HOME)/lib/libhts.a -lz -ldl -lxml2 -lpthread -llzma -lbz2 -lcurl -lcrypto

t_filterfile_SOURCES = t_filterfile.c
t_filterfile_CFLAGS = -DDATA_DIR=$(abs_srcdir)/data

t_posfile_SOURCES = t_posfile.c
t_posfile_CFLAGS = -DDATA_DIR=$(abs_srcdir)/data

t_i2b_SOURCES = t_i2b.c ../src/posfile.c ../src/bclfile.c ../src/filterfile.c ../src/array.c ../src/parse.c ../src/hts_addendum.c
t_i2b_CFLAGS = -I$(top_srcdir)/src -I/usr/include/libxml2 -DDATA_DIR=$(abs_srcdir)/data
t_i2b_LDADD = $(HTSLIB_HOME)/lib/libhts.a -lz -ldl -lxml2 -lpthread -llzma -lbz2 -lcurl -lcrypto

EXTRA_DIST = data

