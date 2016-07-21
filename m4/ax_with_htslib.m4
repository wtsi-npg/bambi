#
# SYNOPSIS
#
#   AX_WITH_HTSLIB()
#
# DESCRIPTION
#
#   This macro searches for the header files and libraries of a htslib
#   (https://github.com/samtools/htslib) installation.  If --with-htslib is specified
#   without argument, or as --with-htslib=yes, a packaged (.deb or
#   .rpm) installation is assumed and the default system paths will be
#   searched.  If --with-htslib=DIR is specified, a run-in-place htslib
#   installation will be searched for in DIR.
#
#   The system header and library paths will be used for run-in-place
#   iRODS installation dependencies, in preference to the those
#   dependencies provided by iRODS in its 'externals' directory
#   because the latter cannot be determined reliably.
#
#   The macro defines the symbol HAVE_HTSLIB if the library is found
#   and the following output variables are set with AC_SUBST:
#
#     HTSLIB_CPPFLAGS
#     HTSLIB_LDFLAGS
#     HTSLIB_LIBS
#
#   You can use them like this in Makefile.am:
#
#      AM_CPPFLAGS       = $(HTSLIB_CPPFLAGS)
#      AM_LDFLAGS        = $(HTSLIB_LDFLAGS)
#      library_la_LIBADD = $(HTSLIB_LIBS)
#      program_LDADD     = $(HTSLIB_LIBS)
#
# LICENSE
#
# Copyright (C) 2016, Genome Research Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

AC_DEFUN([AX_WITH_HTSLIB], [
   HTSLIB_HOME=
   HTSLIB_CPPFLAGS=
   HTSLIB_LDFLAGS=
   HTSLIB_LIBS=

   saved_CPPFLAGS="$CPPFLAGS"
   saved_LDFLAGS="$LDFLAGS"
   saved_LIBS="$LIBS"

   AC_ARG_WITH([htslib],
       [AS_HELP_STRING([--with-htslib[[=DIR]]], [select htslib directory])], [
           AS_IF([test "x$with_htslib" != "xno"], [

               AS_IF([test "x$with_htslib" = "xyes"], [
                   AC_MSG_CHECKING([for packaged htslib])
                   CPPFLAGS="-I/usr/include/htslib $CPPFLAGS"
                   LDFLAGS="-L/usr/lib/htslib/externals $LDFLAGS"
               ], [
                   HTSLIB_HOME="$with_htslib"
                   CPPFLAGS="-I$HTSLIB_HOME/include $CPPFLAGS"
                   LDFLAGS="-L$HTSLIB_HOME/lib $LDFLAGS"
               ])

               HTSLIB_CPPFLAGS="$CPPFLAGS"
               HTSLIB_LDFLAGS="$LDFLAGS"
               HTSLIB_LIBS="$LIBS"

               CPPFLAGS="$saved_CPPFLAGS"
               LDFLAGS="$saved_LDFLAGS"
               LIBS="$saved_LIBS"

               unset saved_CPPFLAGS
               unset saved_LDFLAGS
               unset saved_LIBS

               AC_SUBST([HTSLIB_HOME])
               AC_SUBST([HTSLIB_CPPFLAGS])
               AC_SUBST([HTSLIB_LDFLAGS])
               AC_SUBST([HTSLIB_LIBS])
           ])
   ])
])
