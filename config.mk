#  Optional configure Makefile overrides for bambi
#
#    Copyright (C) 2016 Genome Research Ltd.
#
#    Author: Jennifer Liddle <js10@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This is @configure_input@
#
# If you use configure, this file overrides variables and augments rules
# in the Makefile to reflect your configuration choices.  If you don't run
# configure, the main Makefile contains suitable conservative defaults.


HTSDIR = /software/solexa/pkg/htslib/htslib-1.3
### include $(HTSDIR)/src/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_LIB = $(HTSLIB)
BGZIP = $(HTSDIR)/bin/bgzip
HTSLIB_CPPFLAGS = -I$(HTSDIR)
