/* posfile.c

   Functions to read and parse an Illumina pos, locs, or clocs file.

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
#include <libgen.h>

#include "posfile.h"

posfile_t *posfile_open(char *fname)
{
    posfile_t *posfile = calloc(1, sizeof(posfile_t));
    posfile->current_block = 0;
    posfile->errmsg = NULL;
    posfile->file_type = UNKNOWN_POS;

    // Should we handle compressed (.gz) files?

    char *base = basename(fname);
    char *ext = rindex(base,'.');
    if (ext) ext++;
    if (ext) {
        if (strcmp(ext,"clocs")==0) posfile->file_type = CLOCS;
        if (strcmp(ext,"locs")==0) posfile->file_type = LOCS;
    }
    if (posfile->file_type == UNKNOWN_POS) {
        posfile->errmsg = strdup("posfile_open(): Unknown file type\n");
        return posfile;
    }
    posfile->fhandle = open(fname,O_RDONLY);

    if (posfile->fhandle == -1) {
        posfile->errmsg = strdup(strerror(errno));
        return posfile;
    }

    if (posfile->file_type == CLOCS) {
        read(posfile->fhandle,(void *)&posfile->version,1);
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        read(posfile->fhandle,(void *)&posfile->unread_clusters,1);
        posfile->current_block++;
    }

    if (posfile->file_type == LOCS) {
        // first 8 bytes are unused
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
        read(posfile->fhandle,(void *)&posfile->total_blocks,4);
    }

    return posfile;
}

/*
 * seek to a given cluster number
 */
void posfile_seek(posfile_t *posfile, int cluster)
{
    off_t pos = 12 + cluster * 8;
    if (posfile->file_type != LOCS) {
        fprintf(stderr,"Can only handle NextSeq pos files of type LOC\n");
        exit(1);
    }

    off_t r = lseek(posfile->fhandle, pos, SEEK_SET);
    if (r != pos) {
        fprintf(stderr,"Trying to seek to %d (cluster %d) but returned %d\n", (int)pos, cluster, (int)r);
        perror("posfile_seek() failed");
    }
}

void posfile_close(posfile_t *posfile)
{
    free(posfile->errmsg);
    if (posfile->fhandle >= 0) close(posfile->fhandle);
    free(posfile);
}

static int locs_next(posfile_t *posfile)
{
    float dx, dy;

    if (posfile->current_block >= posfile->total_blocks) return -1;
    posfile->current_block++;

    read(posfile->fhandle, (void *)&dx, 4);
    read(posfile->fhandle, (void *)&dy, 4);

    posfile->x = 10 * dx + 1000.5;
    posfile->y = 10 * dy + 1000.5;

    return 0;
}

static int clocs_next(posfile_t *posfile)
{
    unsigned char dx, dy;

    while (posfile->unread_clusters == 0 && (posfile->current_block < posfile->total_blocks)) {
        if (read(posfile->fhandle, (void *)&posfile->unread_clusters, 1) != 1) return -1;
        posfile->current_block++;
    }

    if (posfile->unread_clusters == 0) return -1;
    posfile->unread_clusters--;

    read(posfile->fhandle, (void *)&dx, 1);    
    read(posfile->fhandle, (void *)&dy, 1);

    posfile->x = 10 * CLOCS_BLOCK_SIZE * ((posfile->current_block - 1) % CLOCS_BLOCKS_PER_LINE) + dx + 1000;
    posfile->y = 10 * CLOCS_BLOCK_SIZE * ((posfile->current_block - 1) / CLOCS_BLOCKS_PER_LINE) + dy + 1000;
    return 0;
}

int posfile_next(posfile_t *posfile)
{
    if (posfile->file_type == CLOCS) return clocs_next(posfile);
    if (posfile->file_type == LOCS) return locs_next(posfile);
    return -1;
}

