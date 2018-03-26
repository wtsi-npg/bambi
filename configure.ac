AC_INIT([bambi], m4_esyscmd_s([git describe --dirty --always --tags]), [js10@sanger.ac.uk])
AC_ARG_VAR(HTSDIR,Directory to look for hts)
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror foreign tar-pax no-dependencies])
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
AC_CHECK_LIB([hts], [sam_hdr_del], [AC_DEFINE([HAVE_SAM_HDR_DEL],[1],[Does htslib contain sam_hdr_del()?])])
PKG_CHECK_MODULES(XML, libxml-2.0 >= 2.7)
CPPFLAGS="$saved_CPPFLAGS"
LDFLAGS="$saved_LDFLAGS"

AC_CHECK_LIB([xml2], [xmlParseFile])
AC_CHECK_LIB([gd], [gdImageCreate])

AC_CONFIG_SRCDIR([src/bambi.h])

AC_CONFIG_FILES([ Makefile ])
AC_OUTPUT

