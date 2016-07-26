AC_INIT([bambi], [1.0], [js10@sanger.ac.uk])
AC_ARG_VAR(HTSDIR,Directory to look for hts)
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_CC

LT_INIT

AC_CONFIG_MACRO_DIR([m4])

AX_WITH_HTSLIB

AC_CONFIG_HEADERS([config.h])

saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
CPPFLAGS="$HTSLIB_CPPFLAGS"
LDFLAGS="$HTSLIB_LDFLAGS"
AC_CHECK_HEADERS([cram/sam_header.h])
AC_CHECK_LIB([hts], [bam_aux_update_str], [AC_DEFINE([HAVE_BAM_AUX_UPDATE_STR],[1],[Does htslib contain bam_aux_update_str()?])])
CPPFLAGS="$saved_CPPFLAGS"
LDFLAGS="$saved_LDFLAGS"

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

