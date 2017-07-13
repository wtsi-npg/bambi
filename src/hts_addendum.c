/* hts_addendum

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bambi.h"
#include "hts_addendum.h"

#ifndef HAVE_BAM_AUX_UPDATE_STR
int bam_aux_update_str(bam1_t *b, const char tag[2], int len, const char *data)
{
    uint8_t *s = bam_aux_get(b,tag);
    if (!s) return -1;
    char type = *s;
    if (type != 'Z') { fprintf(stderr,"bam_aux_update_str() called for type '%c' instead of 'Z'\n", type); abort(); }
    bam_aux_del(b,s);
    s -= 2;
    int l_aux = bam_get_l_aux(b);

    b->l_data += 3 + len;
    if (b->m_data < b->l_data) {
        ptrdiff_t s_offset = (ptrdiff_t) (s - b->data);
        b->m_data = b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
        s = (uint8_t *) (b->data + s_offset);
    }
    memmove(s+3+len, s, l_aux - (s - bam_get_aux(b)));
    s[0] = tag[0];
    s[1] = tag[1];
    s[2] = type;
    memmove(s+3,data,len);
    return 0;
}
#endif

#ifndef HAVE_SAM_HDR_DEL
SAM_hdr *sam_hdr_del(SAM_hdr *hdr, char *type, char *ID_key, char *ID_value) {
    int i,n;
    int *lines;
    char *newtext = malloc(sam_hdr_length(hdr)+15);

    lines = ksplit(&hdr->text,'\n',&n);
    *newtext = 0;
    for (i=0; i<n; i++) {
        char * ln = hdr->text.s + lines[i];
        if (strstr(ln,type) == ln+1) {
            if (!ID_key) continue;
            char *tag = malloc(strlen(ID_key)+2+strlen(ID_value));
            strcpy(tag,ID_key); strcat(tag,":"); strcat(tag,ID_value);
            if (strstr(ln,tag)) { free(tag); continue; }
            free(tag);
        }
        strcat(newtext,ln); strcat(newtext,"\n");
    }
    free(lines);
    sam_hdr_free(hdr);
    hdr = sam_hdr_parse_(newtext,strlen(newtext));
    free(newtext);
    return hdr;
}
#endif

