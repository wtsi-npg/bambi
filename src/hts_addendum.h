/* hts_addendum.h

   Functions which we hope will migrate to htslib one day

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __HTSLIB_ADDENDUM_H__
#define __HTSLIB_ADDENDUM_H__

#include "htslib/sam.h"
#include "cram/sam_header.h"
#include "cram/cram_samtools.h"

#ifndef HAVE_BAM_AUX_UPDATE_STR
int bam_aux_update_str(bam1_t *b, const char tag[2], int len, const char *data);
#endif

#ifndef HAVE_SAM_HDR_DEL
SAM_hdr * sam_hdr_del(SAM_hdr *hdr, char *type, char *ID_key, char *ID_value);
#endif

#endif

