/* bclfile.c

   Functions to read and parse an Illumina BCL file.

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
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "bclfile.h"

#define BCL_BASE_ARRAY "ACGT"
#define BCL_UNKNOWN_BASE 'N'

/*
 * Try to open the given bcl file.
 * If that doesn't work, try appending ".gz" and gzopen it
 * If *that* doesn't work, return the error message in the errmsg field
 */
bclfile_t *bclfile_open(char *fname)
{
    bclfile_t *bclfile = calloc(1, sizeof(bclfile_t));
    bclfile->current_cluster = 0;
    bclfile->gzhandle = NULL;
    bclfile->fhandle = open(fname, O_RDONLY);
    if (bclfile->fhandle == -1) {
        char *gzfname = calloc(1,strlen(fname)+4);
        strcpy(gzfname,fname); strcat(gzfname,".gz");
        bclfile->gzhandle = gzopen(gzfname,"r");
        if (bclfile->gzhandle == NULL) {
            bclfile->errmsg = strdup(strerror(errno));
        } else {
            gzread(bclfile->gzhandle,(void *)&bclfile->total_clusters,4);
        }
    } else {
        read(bclfile->fhandle, (void *)&bclfile->total_clusters, 4);
    }
    return bclfile;
}

void bclfile_close(bclfile_t *bclfile)
{
    if (bclfile->gzhandle) {
        gzclose(bclfile->gzhandle);
    } else {
        close(bclfile->fhandle);
    }
}

int bclfile_next(bclfile_t *bcl)
{
    int c;

    if (bcl->gzhandle) c = gzgetc(bcl->gzhandle);
    else               read(bcl->fhandle, (void *)&c, 1);

    if (c < 0) return c;

    int baseIndex = c & 0x03;   // last two bits
    bcl->quality = (c & 0xfc) >> 2;     // rest of bits
    if (bcl->quality) {
        bcl->base = BCL_BASE_ARRAY[baseIndex];
    } else {
        bcl->base = BCL_UNKNOWN_BASE;
    }

    bcl->current_cluster++;
    return 0;
}

