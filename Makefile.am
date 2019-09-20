ACLOCAL_AMFLAGS= -I m4

AM_CPPFLAGS = -I$(top_srcdir)/src $(HTSLIB_CPPFLAGS) -std=gnu99
AM_CFLAGS = $(XML_CFLAGS) -std=gnu99
AM_LDFLAGS = -rdynamic $(HTSLIB_LDFLAGS)

bin_PROGRAMS = src/bambi src/check_bcl
src_bambi_SOURCES = src/bambi.c \
                    src/bambi.h \
                    src/bambi_utils.c \
                    src/bambi_utils.h \
                    src/decode.c \
                    src/decode.h \
                    src/i2b.c \
                    src/select.c \
                    src/chrsplit.c \
                    src/read2tags.c \
                    src/spatial_filter.c \
                    src/seqchksum.c \
                    src/seqchksum.h \
                    src/adapters.c \
                    src/adapters.h \
                    src/update.c \
                    src/bclfile.c \
                    src/bclfile.h \
                    src/filterfile.c \
                    src/filterfile.h \
                    src/posfile.c \
                    src/posfile.h \
                    src/array.c \
                    src/array.h \
                    src/parse.c \
                    src/parse.h \
                    src/bamit.c \
                    src/bamit.h \
                    src/hash_table.c \
                    src/hash_table.h \
                    src/parse_bam.c \
                    src/parse_bam.h \
                    src/crc.c \
                    src/crc.h

src_check_bcl_SOURCES = src/check_bcl.c src/bclfile.c src/bambi_utils.c src/array.c

nobase_include_HEADERS = src/cram/cram_samtools.h src/cram/pooled_alloc.h src/cram/string_alloc.h

src_bambi_LDADD = $(HTSLIB_LIBS) -lm
src_check_bcl_LDADD = $(HTSLIB_LIBS) -lm

TESTS = test/t_array \
        test/t_bamit \
        test/t_chrsplit \
        test/t_select \
        test/t_bclfile \
        test/t_decode \
        test/t_filterfile \
        test/t_posfile \
        test/t_i2b \
        test/t_read2tags \
        test/t_sf \
        test/t_seqchksum \
        test/t_adapters \
        test/t_update

dist_doc_DATA = README.md LICENSE

check_PROGRAMS = test/t_read2tags \
                 test/t_array \
                 test/t_bamit \
                 test/t_chrsplit \
                 test/t_select \
                 test/t_bclfile \
                 test/t_decode \
                 test/t_filterfile \
                 test/t_posfile \
                 test/t_i2b \
                 test/t_sf \
                 test/t_seqchksum \
                 test/t_adapters \
                 test/t_update

TEST_CFLAGS = -I$(top_srcdir)/src $(XML_CFLAGS) -DDATA_DIR=$(top_srcdir)/test/data
TEST_LDADD = $(HTSLIB_LIBS) -lm

test_t_update_SOURCES = test/t_update.c src/update.c src/bamit.c src/parse.c src/bambi_utils.c src/array.c src/hash_table.c
test_t_update_CFLAGS = $(TEST_CFLAGS)
test_t_update_LDADD = $(TEST_LDADD)

test_t_read2tags_SOURCES = test/t_read2tags.c src/read2tags.c src/array.c src/bamit.c src/parse.c src/bambi_utils.c
test_t_read2tags_CFLAGS = $(TEST_CFLAGS)
test_t_read2tags_LDADD = $(TEST_LDADD)

test_t_select_SOURCES = test/t_select.c src/select.c src/array.c src/bamit.c src/parse.c
test_t_select_CFLAGS = $(TEST_CFLAGS)
test_t_select_LDADD = $(TEST_LDADD)

test_t_chrsplit_SOURCES = test/t_chrsplit.c src/chrsplit.c src/array.c src/bamit.c src/parse.c src/bambi_utils.c
test_t_chrsplit_CFLAGS = $(TEST_CFLAGS)
test_t_chrsplit_LDADD = $(TEST_LDADD)

test_t_bamit_SOURCES = test/t_bamit.c src/bamit.c src/bambi_utils.c
test_t_bamit_CFLAGS = $(TEST_CFLAGS)
test_t_bamit_LDADD = $(TEST_LDADD)

test_t_array_SOURCES = test/t_array.c src/array.c
test_t_array_CFLAGS = $(TEST_CFLAGS)
test_t_array_LDADD = $(TEST_LDADD)

test_t_bclfile_SOURCES = test/t_bclfile.c src/bclfile.c src/array.c src/bambi_utils.c
test_t_bclfile_CFLAGS = $(TEST_CFLAGS)
test_t_bclfile_LDADD = $(TEST_LDADD)

test_t_decode_SOURCES = test/t_decode.c src/array.c src/bamit.c src/hash_table.c src/seqchksum.c src/crc.c src/parse_bam.c src/bambi_utils.c
test_t_decode_CFLAGS = $(TEST_CFLAGS)
test_t_decode_LDADD = $(TEST_LDADD)

test_t_filterfile_SOURCES = test/t_filterfile.c
test_t_filterfile_CFLAGS = $(TEST_CFLAGS)

test_t_posfile_SOURCES = test/t_posfile.c
test_t_posfile_CFLAGS = $(TEST_CFLAGS)

test_t_i2b_SOURCES = test/t_i2b.c src/i2b.c src/posfile.c src/bclfile.c src/filterfile.c src/array.c src/parse.c src/decode.c src/bamit.c src/hash_table.c src/seqchksum.c src/crc.c src/parse_bam.c src/bambi_utils.c
test_t_i2b_CFLAGS = $(TEST_CFLAGS)
test_t_i2b_LDADD = $(TEST_LDADD)

test_t_sf_SOURCES = test/t_sf.c
test_t_sf_CFLAGS = $(TEST_CFLAGS)
#test_t_sf_LDADD = $(TEST_LDADD)

test_t_seqchksum_SOURCES = test/t_seqchksum.c
test_t_seqchksum_CFLAGS = $(TEST_CFLAGS)
#test_t_seqchksum_LDADD = $(TEST_LDADD)

test_t_adapters_SOURCES = test/t_adapters.c src/bamit.c src/bambi_utils.c
test_t_adapters_CFLAGS = $(TEST_CFLAGS)
test_t_adapters_LDADD = $(TEST_LDADD)

EXTRA_DIST = test/data

AM_COLOR_TESTS=always
