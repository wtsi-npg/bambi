AC_INIT([bambi], m4_esyscmd_s([git describe --dirty --always --tags]), [js10@sanger.ac.uk])
AC_CONFIG_MACRO_DIR([m4])
AC_ARG_VAR(HTSDIR,Directory to look for hts)
AM_INIT_AUTOMAKE([subdir-objects foreign tar-pax no-dependencies])
AC_PROG_CC
AM_PROG_CC_C_O
dnl Turn on compiler warnings, if possible
HTS_PROG_CC_WARNINGS
dnl Flags to treat warnings as errors.  These need to be applied to CFLAGS
dnl later as they can interfere with some of the tests (notably AC_SEARCH_LIBS)
HTS_PROG_CC_WERROR(bambi_late_cflags)

LT_INIT

AX_WITH_HTSLIB

AC_CONFIG_HEADERS([config.h])

saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
CPPFLAGS="$HTSLIB_CPPFLAGS"
LDFLAGS="$HTSLIB_LDFLAGS"
AC_CHECK_HEADERS([cram/sam_header.h])
AC_CHECK_LIB([hts], [bam_aux_update_str], [AC_DEFINE([HAVE_BAM_AUX_UPDATE_STR],[1],[Does htslib contain bam_aux_update_str()?])])
AC_CHECK_LIB([hts], [sam_hdr_del], [AC_DEFINE([HAVE_SAM_HDR_DEL],[1],[Does htslib contain sam_hdr_del()?])])
PKG_CHECK_MODULES(XML, libxml-2.0 >= 2.7)
CPPFLAGS="$saved_CPPFLAGS"
LDFLAGS="$saved_LDFLAGS"

AC_CHECK_LIB([xml2], [xmlParseFile])
AC_CHECK_LIB([gd], [gdImageCreate])

AC_CONFIG_SRCDIR([src/bambi.h])

dnl Apply value from HTS_PROG_CC_WERROR (if set)
AS_IF([test "x$bambi_late_cflags" != x],[CFLAGS="$CFLAGS $bambi_late_cflags"])

AC_CONFIG_FILES([ Makefile ])
AC_OUTPUT

