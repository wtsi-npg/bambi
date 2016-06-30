/* bclfile.c

   Functions to read and parse an Illumina BCL file.

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
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <libgen.h>

#include "bclfile.h"

#define BCL_BASE_ARRAY "ACGT"
#define BCL_UNKNOWN_BASE 'N'

/*
 * Try to open the given bcl/scl file.
 * If that doesn't work, try appending ".gz" and gzopen it
 * If *that* doesn't work, return the error message in the errmsg field
 */
bclfile_t *bclfile_open(char *fname)
{
    bclfile_t *bclfile = calloc(1, sizeof(bclfile_t));
    bclfile->current_cluster = 0;
    bclfile->total_clusters = 0;
    bclfile->gzhandle = NULL;
    bclfile->file_type = BCL;
    bclfile->current_base = 0;
    bclfile->filename = strdup(fname);

    // need to find if this is a BCL or SCL file
    char *base = basename(fname);
    char *ext = rindex(base,'.');
    if (ext) ext++;
    if (strcmp(ext,"scl")==0) bclfile->file_type = SCL;
    // FIXME: this isn't going to recognise a .scl.gz file as scl
    // It will probably crash if it doesn't find an extention (ie ext==NULL)

    bclfile->fhandle = open(fname, O_RDONLY);
    if (bclfile->fhandle == -1) {
        char *gzfname = calloc(1,strlen(fname)+8);
        strcpy(gzfname,fname); strcat(gzfname,".gz");
        bclfile->gzhandle = gzopen(gzfname,"r");
        free(gzfname);
        if (bclfile->gzhandle == NULL) {
            strcpy(gzfname,fname); strcat(gzfname,".bgzf");
            bclfile->gzhandle = gzopen(gzfname,"r");
            if (bclfile->gzhandle == NULL) {
                bclfile->errmsg = strdup(strerror(errno));
            }
        } else {
            gzread(bclfile->gzhandle,(void *)&bclfile->total_clusters,4);
        }
    } else {
        read(bclfile->fhandle, (void *)&bclfile->total_clusters, 4);
    }

    return bclfile;
}

void bclfile_seek(bclfile_t *bcl, int cluster)
{
    if (bcl->gzhandle) {
        gzseek(bcl->gzhandle, (z_off_t)(4 + cluster), SEEK_SET);
    } else {
        lseek(bcl->fhandle, (off_t)(4 + cluster), SEEK_SET);
    }
}

void bclfile_close(bclfile_t *bclfile)
{
    if (bclfile->gzhandle) {
        gzclose(bclfile->gzhandle);
    } else {
        close(bclfile->fhandle);
    }
    free(bclfile->filename);
    free(bclfile);
}

int bclfile_next(bclfile_t *bcl)
{
    int i=0;
    static unsigned char c = 0;

    if (bcl->current_base == 0) {
        if (bcl->gzhandle) {
            i = gzgetc(bcl->gzhandle);
            if (i<0) return i;
            c = i;
        } else {
            if (read(bcl->fhandle, (void *)&c, 1) != 1) return -1;
        }
    }

    if (bcl->file_type == SCL) {
        int baseIndex;
        switch (bcl->current_base) {
            case 0: baseIndex = (c >> 6) & 0x03;    break;
            case 1: baseIndex = (c >> 4) & 0x03;    break;
            case 2: baseIndex = (c >> 2) & 0x03;    break;
            case 3: baseIndex = c & 0x03;           break;
        }
        bcl->base = BCL_BASE_ARRAY[baseIndex];
        bcl->current_base++;
        if (bcl->current_base > 3) bcl->current_base = 0;
    } else {
        int baseIndex = c & 0x03;   // last two bits
        bcl->quality = (c & 0xfc) >> 2;     // rest of bits
        if (bcl->quality) {
            bcl->base = BCL_BASE_ARRAY[baseIndex];
        } else {
            bcl->base = BCL_UNKNOWN_BASE;
        }
    }

    if (bcl->current_base == 0) bcl->current_cluster++;
    return 0;
}

