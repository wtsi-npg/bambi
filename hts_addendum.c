/* hts_addendum

   Functions which we hope will migrate to htslib one day

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "hts_addendum.h"

int bam_aux_update_str(bam1_t *b, const char tag[2], int len, uint8_t *data)
{
    uint8_t *s = bam_aux_get(b,tag);
    if (!s) return -1;
    char type = *s;
    if (type != 'Z') { fprintf(stderr,"bam_aux_update_str() called for type '%c' instead of 'Z'\n", type); abort(); }
    bam_aux_del(b,s);
    s -= 2;
    int l_aux = bam_get_l_aux(b);
    uint8_t *aux = bam_get_aux(b);

    b->l_data += 3 + len;
    if (b->m_data < b->l_data) {
        b->m_data = b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
    memmove(s+3+len, s, l_aux - (s - aux));
    s[0] = tag[0];
    s[1] = tag[1];
    s[2] = type;
    memmove(s+3,data,len);
    return 0;
}


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

