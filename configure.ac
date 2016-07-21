AC_INIT([bambi], [1.0], [js10@sanger.ac.uk])
AC_ARG_VAR(HTSDIR,Directory to look for hts)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC

LT_INIT

AC_CONFIG_MACRO_DIR([m4])


AC_SUBST(xCPPFLAGS)

AX_WITH_HTSLIB

AC_CONFIG_HEADERS([config.h])
AC_CHECK_HEADERS([cram/sam_header.h])
AC_CHECK_HEADERS([htslib/sam.h])
AC_SEARCH_LIBS([bam_aux_update_str], [hts])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 test/posfile/Makefile
 test/bclfile/Makefile
 test/filterfile/Makefile
 test/decode/Makefile
 test/i2b/Makefile
])
AC_OUTPUT

