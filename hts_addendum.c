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

#include "hts_addendum.h"

SAM_hdr *sam_hdr_del(SAM_hdr *hdr, char *type, char *ID_key, char *ID_value) {
    int i,n;
    int *lines;
    char *newtext = malloc(sam_hdr_length(hdr)+1);

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

